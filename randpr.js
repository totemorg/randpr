// UNCLASSIFIED

'use strict';
/**
@class RAN
@requires stream
@requires jslab

refs:
www.statslab.cam.ac.uk/~rrw1
www.stat.yale.edu/~pollard
people.math.gatech.edu/~randall
www.stat.berkeley.edu/~pitman
www.math.dartmouth.edu/~pw
**/

var			
	// globals
	TRACE = "R>",

	// nodejs modules
	STREAM = require("stream");			// data streams

const { $, $$, EM, MVN, ME, Copy, Each, Log, FLOW } = require("jslab").libs;
const { sqrt, floor, random, cos, sin, abs, PI, log, exp, min, max} = Math;

class RAN {
	
	constructor(opts, cb) {
		Copy({  // default configuration

			N: 1, 		// ensemble size
			wiener: 0,  // number of steps at each time step to create weiner / SSI process. 0 disables
			symbols: null, 	// state indecies K || ["state", ...] || {state: index, ...} defaults to [0, 1, ...]
			corrMap: null,   // map state index to correlation value [value, ... ]
			jumpModel: "",   // inhomogenous model (e.g. "gillespie" ) or "" for homogeneous model 
			store: 	null,  // created by pipe()
			steps: 1, // number of process steps of size dt 
			ctmode: false, 	// true=continuous, false=discrete time mode 
			obslist: null, // observation save list
			keys: null,  // event key names
			
			learn: null, 	// sup/unsup learner(supercb) with callbacks supercb(evs) || supercb(null,onend)
						
			// K-state config parameters
			trP: null, 	// [KxK] from-to state trans probs or { states:K, index: { index: prob, ...}, ...}
			emP: null, // {mu: mean, sigma: stddevs, dims: [dim, ....] } xyz-emmision probs

			// K=2 state config parameters
			alpha: 0,  // on-to-off rate [jumps/s]
			beta: 0,  // off-to-on rate [jumps/s]
			p: 0,  // on state probability
			q: 0,  // off(not on) state probability

			filter: function (str,ev) {  // filter output event ev to store/stream str
			/**
			Output event filter
				filter: function (str, ev) { // event ev for stream/store str
						switch ( ev.at ) {   // streaming plugins provide an "at" to filter events on
							case "...":
							case "...":
								str.push(ev);	// return the event
						}
					}  
			*/
				str.push( ev ); 
			},

			// internal variables
			K: 0, 		// number of states
			U: null,    // [N] states at time t
			U0: null, 	// [N] initial states at time t = 0
			U1: null, 	// [N] state buffer 
			J: null,  // [N] jump counter
			HT: null, 	// [N] next jump time
			WU: null, 		// [N] wiener/SFI ensemble
			WQ: null, 		// [N] wiener/SFI ensemble for cummulative walks
			
			RT: null, 	// [KxK] from-to holding (mean recurrence) times
			abT: null, 	// [K'] absorption times K' <= K
			abP: null,	// [K' x K-K'] absorption probabilities K' <= K
			mleA: null, 	// [KxK] from-to state mle trans probabilities
			mleB: null, 	// {mu,sigma} observation mixing parameters
			corP: null, 	// [KxK] stat correlation probabilities
			cumP: null,	// [KxK] from-to cummulative state transition probabilities
			N0: null, 	// [KxK] from-to cummulative counts in to-state given starting from-state
			N1: null,	// [KxK] from-to state transition counts
			cumH: null, 	// [KxK] cummulative time in from-to transition
			cumN: null, 	// [KxK] cummulative number of from-to jumps
			eqP: null, 	// [K] equilibrium state probabilities 
			A: null,	// [KxK] jump rates
			
			// supervised learning parms			
			batch: 0, 				// batch size in dt-steps to make MLEs
			
			// sampling parms
			halt: false, // default state when learning
			Tc: 0,  // coherence time >0 [dt] 
			t: 0, 	// time [dt]
			s: 0, 	// step count
			dt: 1, 	// sample time 
			jumps: null, // [maxJumps] distribution of state jumps for event counting
			samples: 0 // number of ensemble members sampled
		}, this);

		if (opts) Copy(opts, this);
		
		var 
			ran = this,
			N = this.N, // ensemble size
			trP = this.trP, // transition probs KxK or coersed to KxK
			emP = this.emP, // emission (aka observation) probs
			keys = this.keys = Copy(this.keys || {}, { index:"n", state:"u", class:"k", x:"x", y:"y", z:"z", t:"t" }), // event keys
			//obs = this.obs, // observation (aka emission or mixing) parms
			symbols = this.symbols || 2;  // state symbols

		if ( this.p ) {   // two-state process via p,Tc
			var
				p = this.p,
				Tc = this.Tc, 
				meanrate = 2.3/Tc,
				rate = 2*meanrate,
				alpha = p * rate,
				beta = rate - alpha,
				q = beta / rate,
				K = this.K = 2;
			
			trP = this.trP = [[1-p, p], [q, 1-q]];
		}

		else
		if ( this.alpha )  { // two-state process via alpha,beta
			var 
				alpha = this.alpha,
				beta = this.beta, 
				p = alpha / (alpha + beta), 
				q = beta / (alpha + beta),
				K = this.K = 2;
			
			trP = this.trP = [[1-p, p], [q, 1-q]];
		}

		if ( emP ) {
			this.obslist = [];
			if (dims = emP.dims) {
				var 
					K = 1,
					drop = dims.use( (n,Dims) => K *= Dims[n] ),
					weights = emP.weights,
					D = dims.length,
					grid = emP.grid = perms( [], dims, []),  // state grid	
					mus = emP.mu = [],
					sigmas = emP.sigma = [],
					gen = emP.gen = $(K, function (k,Gen) { // gauss mixing (mu,sigma) parms 
						var 
							n = 0,

							mu = $(D, (i,mu) =>
								mu[i] = grid[k][n++] + 0.5
							),

							L = $$(D,D, (i,j, L) => 	// lower trianular matrixfs with real, positive diagonal
								L[i][j] = (i <= j ) ? random() : 0
							), 

							sigma = $$(D,D, function (i,j, A) { // hermitian pos-def matrix via cholesky decomp
								var dot = 0;
								L.use( function (n) {
									dot += L[i][n] * L[j][n];
								});
								A[i][j] = dot * weights[i] * weights[j]
							});

						mus.push( mu );
						sigmas.push( sigma );
						
						Gen[k] = MVN( mu, sigma ).sample;
					});
				
				this.K = K;
			}
			
			else {
				var
					mus = emP.mu,
					sigmas = emP.sigma,
					K = this.K = mus.length,
					gen = emP.gen = $(K, (k,Gen) => Gen[k] = MVN( mus[k], sigmas[k] ).sample );
			}
		}

		//if (!symbols) symbols = K;
		
		switch ( symbols.constructor.name ) {
			case "Object":
				K = 0;
				for (var key in symbols) K++;
				this.K = K;
				break;
				
			case "Array":
				for (var syms = {}, k=0, K=symbols.length; k<K; k++) syms[ symbols[k] ] = k;
				this.K = K;
				symbols = this.symbols = syms;
				break;
				
			default:
				for (var syms = {}, k=0, K=symbols; k<K; k++) syms[ k ] = k;
				this.K = K;
				symbols = this.symbols = syms;
		}

		switch ( trP.constructor.name ) {
			case "Function":
				break;
				
			case "String":
				var K = this.K = symbols.length;
				switch (trP) {
					case "random":
						trP = this.trP = $$(K, K, (fr,to,P) => P[fr][to] = random() );
						trP.use( function (fr) {
							var 
								P = trP[fr],
								sum = P.sum( );

							P.use( (to,P) => P[to] /= sum);
						});
						break;
				}
				
				break;

			case "Array":
				var
					K = this.K = trP.length;

				break;
				//Kto = trP[0].length,
				//P = $$(K, K, (fr,to,P) => P[fr][to] = (fr<Kfr && to<Kto) ? trP[fr][to] : 0 );
				//trP = this.trP = P;

			case "Object":
				var
					K = this.K = trP.states || 2,
					P = $$(K, K, $$zero),
					dims = emP ? emP.dims : [K];

				delete trP.states;
				for (var frKey in trP) {
					var 
						frP = trP[frKey],
						frIndex = index( frKey.split(","), dims );

					//Log("fr", frKey, frIndex);

					for (var toKey in frP) {
						var toIndex = index( toKey.split(","), dims );
						P[frIndex][toIndex] = frP[toKey];
					}
				}

				balanceProbs(P);
				trP = this.trP = P;
				//Log("trP", trP);
				break;
		}
		
		var 
			RT = this.RT = meanRecurTimes(trP),  // from-to mean recurrence times
			ab = this.ab = firstAbsorb(trP),  // first absoption times, probs, and states
			eqP = this.eqP = $(K, (k,eqP) => eqP[k] = 1/RT[k][k]	);  // equlib state probs
		
		if ( trP.constructor == Array ) 
			if ( trP[0].constructor == Array ) 
				this.transMode = "homogen";
			else
				this.transMode = "mh";
		else
			this.transMode = "external";
		
		//Log(K, trP, N);
		
		// create zero-mean correlation map
		var map = this.corrMap = new Array(K);
		
		if (K % 2) {
			map[0] = 0;
			for (var a=1, k=1; k<K; a++) {
				map[k++] = a; 
				map[k++] = -a;
			}
		}

		else			
			for (var a=1, k=0; k<K; a++) {
				map[k++] = a; 
				map[k++] = -a;
			}

		Log(TRACE, {
			keys: keys,
			states: K, 
			syms: symbols, 
			xMap: map,
			trProbs: trP
		});	

		// allocate the ensemble
		var 
			U1 = this.U1 = $(N),
			J = this.J = $(N, $zero),
			N1 = this.N1 = $$(K,K,$$zero),	
			mleA = this.mleA = $$(K,K,$$zero), 
			cumH = this.cumH = $$(K,K,$$zero),
			cumN = this.cumN = $$(K,K,$$zero),
			cumP = this.cumP = $$(K,K,(fr,to,P) => P[fr][to] = trP[fr][to] ),
			Rmle = this.Rmle = $$(K,K),
			err = this.err = 1,
			corP = this.corP = $$(K,K),
			//obs=this.obs,
			emP = this.emP,
			p = 1/K,
			Np = p * N,
			N0 = this.N0 = $$(K,K, (fr,to,N0) => N0[fr][to] = delta(fr,to) ? Np : 0 ),
			/*
			this.learn 
				? $$(K,K,$$zero)
				: $$(K,K, (fr,to,N0) => N0[fr][to] = delta(fr,to) ? Np : 0 ),
			*/
			
			HT = this.HT = $(N),
			U = this.U = $(N),
			U0 = this.U0 = $(N),
			ctmode = this.ctmode,
			WU = this.WU = this.wiener ? $(N) : [],
			WQ = this.WQ = this.wiener ? $(N) : [];
		
		this.t = this.s = this.samples = 0;  // initialize process counters

		// initialize ensemble
		
		if ( false ) { // this.learn ) {  // in learning mode
			U.use( (n) => HT[n] = U0[n] = U[n] = 0 );
		}
		
		else { // generative mode
			cumP.use(  (fr) => {  // initialize cummulative state transition probabilties
				cumulative( cumP[fr] );  
			});

			if (K == 2) {  // initialize two-state process (same as K-state init but added control)
				var R01=RT[0][1], R10=RT[1][0];

				U.use( (n) => {
					if ( n < Np ) {
						var fr = U0[n] = U[n] = 1;
						HT[n] = RT[fr][fr] = ctmode ? expdev(-1/A[fr][fr]) : 0;
					}

					else {
						var fr = U0[n] = U[n] = 0;
						HT[n] = RT[fr][fr] = ctmode ? expdev(-1/A[fr][fr]) : 0;
					}
				});
			}

			else  // initialize K-state process
				U.use( (n) => {
					var fr = floor(random() * K);
					U0[ n ] = U[n] = fr; 
					HT[ n ] = RT[fr][fr] = ctmode ? expdev(-1/A[fr][fr]) : 0;
				}); 
		}
		
		//Log("HT", HT);
		if ( this.wiener ) {  //  initialilze wiener processes
			this.NRV = MVN( [0], [[1]] );
			for (var n=0; n<N; n++) WU[n] = WQ[n] = 0;
		}

		this.gamma = $(this.steps, $zero);
		
		if (cb) cb(null);
	}
	
	statCorr( ) {  // statistical correlation function
		var 
			K = this.K, map = this.corrMap, cor = 0, corP = this.corP, p, N0 = this.N0, N = this.N, samples = this.samples;

		if (samples)
			map.use( (fr) => {
				map.use( (to) => {
					p = corP[fr][to] = N0[fr][to] / samples;
					cor += map[fr] * map[to] * p;
				});
			});
		
		else
			cor = 1;
		
		this.samples += N;

		return cor ; 
	}

	step (evs, cb) {  // advance process forward one step (with events evs if in learning mode)
		
		function draw( P ) { // draw random state with cumulative prob P
			var to = 0, K = P.length;

			for (var u = random(); to < K && P[to] <= u; to++) ;
			return (to == K) ? to-1 : to;
		}

		var 
			ran = this,
			U=this.U,HT=this.HT,RT=this.RT,U0=this.U0,mleA=this.mleA,N0=this.N0,K=this.K,
			cumP = this.cumP, eqP = this.eqP, trP = this.trP,
			U1=this.U1, N1=this.N1, cumH = this.cumH, cumN = this.cumN, A = this.A, 
			symbols=this.symbols, keys = this.keys, emP = this.emP, J=this.J,
			t = this.t, N = this.N, s=this.s,
			trans = {
				external: trP,
				
				homogen: function ( fr ) {  // homogeneous transitions
					var to = draw( cumP[fr] );
					return to;
				},

				// inhomogeneous jumps

				gillespie: function( fr ) {  // return toState by computing cumulative trans probs P based on holding times RT
					var 
						P = cumP[fr],
						R0 = RT[fr], 
						K = P.length;

					P.use( (to) => P[to] = (fr==to) ? 0 : RT[to] / R0 );

					cumulative(P);	
					var P0 = P[K-1];
					P.use( (to) => P[to] /= P0 );

					return draw( P[fr] );
				},

				mcmc: function ( fr ) {
				},

				mh: function ( fr ) {  // return toState using metropolis-hastings with generator G
					var 
						P = eqP,
						toG = cumP[fr], 
						to = draw(toG),
						frG = cumP[to],
						Ap = P[to] / P[fr],
						Ag = frG[to] / toG[fr],
						accept = min(1 , Ap * Ag), // acceptance criteria
						u = random();

					return (u <= accept) ?  to : fr;
				}
			},
			stateTrans = trans[ran.transMode];
				
		U.use( (n) => U1[n] = U[n] );  // hold states to update the N1 counters
		
		this.gamma[s] = this.statCorr();
		
		if (evs) { // in learning mode with time-ordered events
			/*
			if ( !t ) { // initialize states at t=0
				ran.gamma[0] = 1;  // force
				evs.forEach( (ev) => {
					var n = ev[keys.index];		// ensemble index
					U[ n ] = symbols[ev[keys.state]] || 0;  // state (0 if hidden)
				});
				U.use( (n) => {
					U0[n] = U1[n] = U[n];		// initialize states
					J[n] = -1;  // force initial jump counter to 0
				}); 
			} */
				
			t = this.t = evs[0].t;	
			//Log("t",t);
			evs.forEach( (ev) => {   // assume events are time-ordered 
				var 
					n = ev[keys.index],  // ensemble index = unique event id
					fr = U[n],   // current state of member (0 if hidden)
					to = symbols[ev[keys.state]];	// event state (if supervised) or undefined (if hidden)
				
				cumH[fr][to] += t - HT[n]; 	// total holding time in from-to jump
				cumN[fr][to] ++;  	// total number of from-to jumps

				J[ n ]++; 		// increment jump counter
				U[ n ] = to || 0;  		// update state (0 if hidden)
				HT[ n ] = t;			// hold current time
				
				ran.onEvent(n,to,0, [ev[keys.x], ev[keys.y], ev[keys.z]] ); 	// callback with jump info
			});
		}

		else   // generative mode
			U.use( (n) => {
				
				var
					frState = U[n],
					toState = stateTrans( frState );
				
				if ( frState != toState) { // jump if state changed
					var
						held = t - HT[n],	// initially 0 and remains 0 in discrete-time mode
						hold = this.ctmode ? expdev( 1/A[frState][toState] ) : 0 ;  // draw expected holding time
					
					cumH[frState][toState] += held; // cummulative holding time in from-to jump
					cumN[frState][toState] ++;  // cummulative number of from-to jumps
					RT[frState][frState] = hold;  // update expected holding time 
					
					J[ n ]++; 		// increment jump counter
					U[ n ] = toState;  		// set new state
					HT[ n ] = t + hold;    // advance to next jump time (hold is 0 in discrete time mode)
						
					ran.onEvent(n, toState, hold, emP ? emP.gen[toState]() : null ); 	// callback with jump info
				}
				
			});	
			
		//Log( (t<10) ? "0"+t : t, U.join(""));
		//if (t<50) Log( t<10 ? "0"+t : t,U,J);		
		//if (t<5) Log(t,N0);
		
		U.use( (n) => {   // adjust initial from-to counters for computing ensemble correlations
			N0[ U0[n] ][ U[n] ]++; 
		});

		U.use( (n) => {  // adjust step from-to counters for computing trans probs
			N1[ U1[n] ][ U[n] ]++;
		});

		if ( ran.wiener ) {  // step wiener process
			//Log("wsteps", this.wiener);
			var 
				M = ran.wiener,
				//floor = Math.floor, sqrt = Math.sqrt, 
				nrv = ran.NRV.sample, WU = ran.WU, WQ = ran.WQ;
			
			for (var n=0; n<N; n++) {				
				for (var sum=WQ[n], j=1, J=floor(M*t); j<=J; j++) sum += nrv()[0];
				WU[n] = sum / sqrt(M);
				WQ[n] = sum;
			}
		}
		
		ran.onStep();
		ran.t += ran.dt;
		ran.s++;
	}
	
	start ( ) {	  // start process in learning (reverse) or generative (forward) mode
		var 
			ran = this,
			U = this.U,
			trP = this.trP,
			batch = this.batch;

		if ( ran.learn && !ran.halt )  // learning mode
			ran.learn( function supervisor(evs, cb) {  // process events when evs, or terminate with callback(results) when evs exhausted

				if (evs) {
					//Log("FEEDING "+evs.length + " len="+evs[0].t);
					ran.step(evs);
				}

				else {
					//Log("HALTING", ran.t, ran.steps);
					ran.halt = true;
					ran.onEnd();
					if (cb)
						cb({  // callback with a ran ctx 
							trP: ran.trP,	// transition probs
							store: ran.store,  // output event store
							//steps: ran.steps,	// process steps
							T: ran.steps,		// observation time
							F: ran.F,	// event count frequencies
							J: ran.J,		// ensemble jump counts
							N: ran.N		// ensemble size
						});
				}

				if ( batch )
					if ( ran.s % batch == 1 ) ran.onBatch();
			});
		
		else { // generative mode
			//Log("start gen", ran.steps, ran.N);
			//U.use( (n) => ran.onEvent (n,U[n],0,Y[n]) );
			
			//ran.step();
			//if ( batch ) ran.onBatch();
			
			while (ran.s < ran.steps) {  // advance process to end
				ran.step(null);
				
				if ( batch )
					if ( ran.s % batch == 1 ) ran.onBatch();
			}
			
			ran.onEnd();
		}
		
	}
	
	corrTime ( ) {  // return correlation time computed as area under normalized auto correlation function
		for (var t=0, T = this.t, Tc = 0; t<T; t++) Tc += abs(this.gamma[t]) * (1 - t/T);
		
		return this.Tc = Tc * this.dt / this.gamma[0] / 2;
		Log(">>>>>>>>>Tc=", Tc);
	}
	
	record (at, ev) {  // record event ev labeled at to store or stream
		ev.t = this.t;
		ev.at = at;
		this.filter(this.store, ev);
	}
	
	onBatch () {    // record MLE jump rates and trans probs
		var 
			ran = this,
			K = this.K,
			Kbar = this.J.avg(),
			t = this.t,
			s = this.s,
			cumH = this.cumH, 
			cumN = this.cumN,
			RT = this.RT,
			trP = this.trP,
			Rmle = this.Rmle,
			N1 = this.N1,
			max = Math.max,
			obslist = this.obslist,
			mleA = this.mleA;
		
		Rmle.use( (fr,to) => {   // estimate jump rates using cummulative HT[fr][to] and N[fr][to] jump times and counts
			Rmle[fr][to] = delta(fr,to) ? 0 : cumH[fr][to] / cumN[fr][to];
		});
		
		N1.use( (fr) => {  // estimate transition probs using the 1-step state transition counts
			var Nfr = N1[fr], Afr = mleA[fr];
			Nfr.sum( (sum) => {
				Afr.use( (to) => {
					Afr[to] = Nfr[to] / sum;
					//Log(fr,to,Nfr[to], sum);
				});
			});
		});
		
		var 
			n = 0,
			err = this.err = trP   // relative error between mle and actual trans probs
				? abs( mleA[n][n] - trP[n][n] ) / trP[n][n]
				: 0;
			
		/*
		$$use(Rmle, function (fr,to) {
			err[fr][to] = (fr == to) ? 0 : ( Rmle[fr][to] - RT[fr][to] ) / RT[fr][to] ;
		}); */
		
		//Log( T, mleA);
		this.record("batch", {
			rel_error: err,
			mle_em_events: obslist ? obslist.length : 0,
			mle_tr_probs: mleA,
			stat_corr: this.gamma[ s-1 ]
		});		
	}

	onError( msg ) {	// record process error condition
		Log(msg);
		this.record("error", { 
			error: msg
		});
	}
	
	onEvent (index,state,hold,obs) {  // record process event info
		
		var obslist = this.obslist;
		
		if (obslist) obslist.push( obs );  // retain for training Viterbi emission probs
		
		this.record("jump", {
			index: index, state:state, hold:hold, obs:obs
		});
	}

	onStep () {		// record process step info
		this.record("step", {
			gamma:this.gamma[this.s],
			walk: this.wiener ? this.WU : []
		});
	}

	onConfig() {  // record process config info
		var
			//obs = this.obs,		
			trP = this.trP,
			emP = this.emP,
			K = this.K,
			RT = this.RT, // = meanRecurTimes(trP),  // from-to mean recurrence times
			eqP = this.eqP = $(K, (k,eqP) => eqP[k] = 1/RT[k][k] ),   // eq state probs
			ab = this.ab = firstAbsorb(trP);			
		
		this.record("config", {
			states: this.K,
			ensemble_size: this.N,		
			sample_time: this.dt,
			nyquist: 1/this.dt,
			cum_tr_probs: this.cumP,
			tr_probs: trP,
			mean_recurrence_times: RT,
			eq_probs: eqP,
			initial_activity: this.p,
			wiener_walks: this.wiener ? "yes" : "no",
			mixing: emP,
			/*mixing: obs 
				? {
					mu: obs.mu,
					sigma: obs.sigma
				} 
				: null, */
			//jump_rates: this.A,
			//avg_jump_rate: this.lambda,
			//exp_coherence_time: this.Tc,
			run_steps: this.steps,
			absorb: ab
		});
	}
	
	onEnd() {  // record process termination info
		
		//Log("onend", this.obslist.length);
		
		var 
			ran = this,
			batch = ran.batch,
			T = ran.steps,
			Tc = ran.corrTime(),
			Kbar = ran.J.avg(),
			M = T / Tc,
			delta = Kbar / M,
			J = ran.J,  // number of state jumps (aka events or counts) made by each ensemble member
			F = ran.F = $( J.max()+1, $zero ),  // count frequencies
			obslist = this.obslist,		
			K = this.K,
			mleB = this.mleB = obslist ? EM( obslist, K) : null;

		J.use( (n) =>F[ J[n] ]++ ); // count frequencies across the ensemble

		//Log("onend J", J);
		
		this.record("end", {  // record supervised stats
			stats: {
				mle_holding_times: ran.Rmle,
				rel_error: ran.err,
				mle_em_probs: ran.mleB,
				mle_tr_probs: ran.mleA,
				tr_counts: ran.N1,
				mean_count: Kbar, 
				coherence_time: Tc, 
				coherence_intervals: M,
				correlation_0lag: ran.gamma[0],
				mean_intensity: Kbar / T,
				degeneracy_param: delta,
				snr: sqrt( Kbar / (1 + delta ) )
			}
		});

		if (evs = this.store)
			evs.forEach( (ev) => {
				ev.s = ev.t / Tc;
			});
	}

	end(stats, saveStore) {  // terminate process
		this.record("end", {  // post learning stats
			stats: stats ? Copy(stats,{}) : {error:"stats unavailable"} 
		});  
		if (saveStore) saveStore( this.store );
	}
	
	pipe(sinkStream) {  // pipe events to a sinking stream or to a callback sinkStream(events)
		var 
			ran = this,
			sync = (typeof sinkStream) == "function";

		Trace( `PIPE${sync ? "sync" : "async"}` );
		
		ran.store = sync
			? []
			: new STREAM.Readable({  // prime and terminate the pipe
				objectMode: true,
				read: function () {  // prime or terminate the pipe
					//Log("randpr pipe at", ran.t);

					if ( ran.s < ran.steps ) 	// prime
						ran.start( );

					else  { // terminate
						ran.end();
						this.push(null);
					}
				}
			});
		
		ran.onConfig();		// process configured
		
		if  (sync) {  // pipe is sync mode using array store
			ran.start();
			
			sinkStream( ran.store );
		}
			
		else {	// pipe in async mode
			var
				ranStream = ran.store,

				editStream = new STREAM.Transform({  // 2nd stage filters events
					writableObjectMode: true,
					readableObjectMode: true,
					transform: function (ev,en,cb) {
						ran.filter(this, ev, this.learn);
						cb(null);
					}
				}),

				charStream = new STREAM.Transform({  // 3rd stage makes events human readable 
					writableObjectMode: true,
					transform: function (ev,en,cb) {
						this.push( JSON.stringify(ev) ); 
						cb(null);
					}
				});

			ranStream.pipe(editStream).pipe(charStream).pipe(sinkStream);
		}
			
	}
		
}

module.exports = RAN;

function expdev(mean) {
	return -mean * log(random());
}

function avgRate(A) {  // computes average jump rate in A not necessarily balanced

	for (var fr=0,lambda=0,K=A.length; fr<K; fr++)
		for (var to=0,Afr=A[fr]; to<K; to++)
			if ( fr != to ) lambda += Afr[to];

	return lambda / (K*K-K); 	
}

function $$zero(i,j,A) {
	A[i][j] = 0;
}

function $zero(i,A) {
	A[i] = 0;
}

function cumulative( P ) {  // replace P with its cumulative
	switch (0) {
		case 0:
			P.use( (k) => {
				if (k) P[k] += P[k-1];
			});
			break;
			
		case 1:
			var cum = 0;
			P.use( (k) => {
				var hold = P[k];
				P[k] = cum;
				cum += hold;
			});
			break;
	}
}

function range (min,max) { // unused - generate a range
	var rtn = new Array(max-min+1);
	for (var n=min,m=0,M=rtn.length; m<=M; m++) rtn[m] = n += 1;
	return rtn;
}	

function balanceRates(A) {   // enforce global balance on jump rates
	A.use( (k) => A[k][k] = - A[k].sum() );
	return A;
}

function balanceProbs(P) {  // enforce global balance on probs
	P.use( (k) => {
		P[k][k] = 1 - P[k].sum() 
	});
	return P;
}			

function delta(fr,to) {
	return (fr==to) ? 1 : 0;
}

function Trace(msg) {
	TRACE.trace(msg);
}

function firstAbsorb(P) {  //< compute first absorption times, probs and states
	var 
		K = P.length,
		kAb = [],
		kTr = [],
		x = P.use( (k) => {
			if ( P[k][k] == 1 ) 
				kAb.push(k+1);
			else
				kTr.push(k+1);
		}),
		ctx = {
			P: ME.matrix(P),
			K: K,
			kAb: ME.matrix(kAb),
			kTr: ME.matrix(kTr),
			nAb: kAb.length,
			nTr: kTr.length,
			abT: ME.matrix([]),
			abP: ME.matrix([])
		};
	
	//Log("ab ctx", JSON.stringify(ctx));
	if ( ctx.nAb && ctx.nTr )
		ME.eval("Q = P[kTr,kTr]; RT = P[kTr,kAb]; N = inv( eye(nTr,nTr) - Q ); abT = N*ones(nTr,1); abP = N*RT;", ctx);
		
	return {
		times: ctx.abT._data,
		probs: ctx.abP._data,
		states: kAb
	};
}

function meanRecurTimes(P) {  //< compute mean recurrence times
/*
If the process is/were Regular, we could itterate the process (e.g. compute some power of the from-to 1-step transition $$ P) to determine 
the equlibrium 
prob vector w, and therefore its associated eq prob $$ W = [w ; w; ... ].  However, in general, P is not Regular.  We require, however, 
that the process P be at 
least Ergodic (w or w/o absorbing states) and, thus, it must possess mean recurrence times H.  So while the computed H must have nontrivial
values for an absorbing P, there is (of course, and by definition) no guarantee that all states will be hit, and thus there	is no guarantee that 
the MLE H will match the computed H at transitions that are never hit.  So, in the general ergodic case, the equib probs w must be determined
by examining the left-nullspace for ( I - P ) whose inverse does not exists (see [*]).  There is, however, an alternative way to compute the w since
sum(w) = 1 by definition.  Thus we decompose ( I - P ) as follows:

		w * P = w
		[ 1, wk ] [ [P0 Pu] ; [Pl Pk] ] = [1, wk]

leaving 2 simultaneous equations:

		P0 + wk * Pl = 1
		Pu + wk * Pk = wk

the last of which is solved for a wk, thence w = renormalized( [1,wk ] ) such that sum(w) = 1.

This technique fails, however, when det(Pk - I ) vanishes, that is, when wk (and therefore w) is not 
unique; we test for 
such non-ergodic P by testing the det( Pk - I ). [*] Note further that in the higher KxK space, det(P - I) 
always
vanishes as it must (the columns of P must sum to 1).  Indeed, if (P - I) has an inverse A, then 
P = I + inv(A);
but inv(A) cannot exists as A is balanced (its columms summing to 0).  Thus, P does not uniquely 
determine the process: only the mean recurrence times H and the equlib pr w determine the process.
*/	

	var 
		ctx = {
			P: ME.matrix(P),
			K: P.length
		},
		K = ctx.K;

	if ( K > 1) {
		ME.eval("k=2:K; P0=P[1,1]; Pl=P[k,1]; Pu=P[1,k]; Pk=P[k,k]; A = Pk - eye(K-1); Adet = abs(det(A)); ", ctx);

		Log(TRACE, {"MRT det": ctx.Adet});

		if ( ctx.Adet < 1e-3 ) {
			Log("Proposed process is not ergodic, thus no unique eq prob exist.  Specify one of the following eq state prs: P^inf --> ", ME.pow(P,20));
			return $$(K,K, $$zero );
		}

		else {
			ME.eval("wk= -Pu*inv(A);", ctx);

			ctx.w = ME.matrix([1].concat(ctx.wk._data[0]));

			ME.eval("w = w / sum(w); w = [w]; Z = inv( eye(K) - P + w[ ones(K) , 1:K] ); H = zeros(K,K); ", ctx);

			var 
				H = ctx.H._data,
				Z = ctx.Z._data,
				w = ctx.w._data[0];

			for (var fr=0;fr<K; fr++) 
				for (var to=0; to<K; to++) 
					H[ fr ][ to ] = ( ( fr == to ) ? 1 / w[ to ] : ( Z[ to ][ to ] - Z[ fr ][ to ] ) / w[ to ] );	

			return H;
		}
	}
	
	else
		return [[1]];
}

function perms(vec,dims,vecs,norm) {  //< generate permutations

	if (vec.length == dims.length) 
		vecs.push(vec);
	
	else 
		for (var idx = 0, max = dims[vec.length]; idx<max; idx++) 
			perms(vec.concat( norm ? norm(idx,max) : idx), dims, vecs,norm);
	
	return vecs;
}

function poisson(m,a) {
	// a^m e(-a) / m!
	for (var sum=0,k=m; k; k--) sum += log(k);
	return exp( m*log(a) - a - sum );	
}

function dirichlet(alpha,grid,logP) {  // dirchlet allocation
	var 
		K = alpha.length,
		N = x[0].length,
		logBs = $(K, (k,B) => B[k] = GAMMA.log( alpha[k] ) ),
		logB = logBs.sum() - GAMMA.log( alpha.sum() );
	
	grid.use( (n) => {
		var
			logAs = $(K, (k,A) => A[k] = (alpha[k] - 1) * log( grid[k] ) ),
			logA = logAs.sum();
	
		logP[n] = logA - logB;
	});
}	

function index(keys, dims) {
	var N  = dims.length, idx = 0, off = 1;
	
	if (keys.length == 1)
		idx = parseInt(keys[0]);
	
	else
		for (var n=0; n<N; off *= dims[n], n++) 
			idx += off * parseInt( keys[n] );

	//Log( keys, idx);
	
	return idx;
}

[ 
	function sample(delta) {
		var 
			A = this,
			k = 0,
			rtn = $( floor(A.length/delta), (n,RT) => RT[n] = A[k += delta] );
			return rtn;
	},
	
	function sum(cb) {
		for (var A=this, Sum=0, k=0, K= A.length; k<K; k++) Sum+= A[k];

		if (cb) cb(Sum,this);

		return Sum;
	},

	function avg() {
		return this.sum() / this.length;
	},

	function max() {
		var A = this, Amax = -1e99, Aidx = 0;
		A.use( (k) => {
			if ( A[k] > Amax ) {
				Amax = A[k];
				Aidx = k;
			}
		});
		return Amax;
	}

	/*function use(cb) {	// use vector A with callback cb(idx,A)
		var A = this, N = A.length;

		if (A.rows) {
			var M = A.rows, N = A.columns;

			for (var m=0; m<M; m++) for (var n=0, Am = A[m]; n<N; n++) cb(m,n,A,Am);
			return A;
		}

		else
			for (var n=0,N=A.length; n<N; n++) cb(n,A);

		return A;
	}	*/
].extend(Array);

switch (process.argv[2]) {   //< unit tests
	case "R1":  // mean recurrence times
		Log( meanRecurTimes(  
			[[0.5,0.25,0.25],[0.5,0,0.5],[0.25,0.25,0.5]]   // regular and ergodic
/*
MRT det= 0.375
[ [ 2.5, 4, 3.3333333333333335 ],
  [ 2.6666666666666665, 5, 2.6666666666666665 ],
  [ 3.3333333333333335, 4, 2.5 ] ]

*/
			// [[0,1,0,0,0], [.25,0,.75,0,0], [0,.5,0,.5,0], [0,0,.75,0,.25], [0,0,0,1,0]]  // not regular but ergodic (no absorbing states)
/*
MRT det= 0.09375
[ [ 16.000000000000007,
1,
2.666666666666666,
6.333333333333334,
21.333333333333336 ],
[ 15.000000000000004,
4,
1.6666666666666665,
5.333333333333334,
20.333333333333336 ],
[ 18.66666666666667,
3.666666666666667,
2.666666666666666,
3.666666666666667,
18.666666666666668 ],
[ 20.33333333333334, 5.333333333333334, 1.666666666666666, 4, 15 ],
[ 21.33333333333334, 6.333333333333334, 2.666666666666666, 1, 16 ] ]

*/
		)); 
		break;
		
	case "R2":	  // absorption times
		Log( firstAbsorb( 
			//[[1,0,0,0,0],[0.5,0,0.5,0,0],[0,0.5,0,0.5,0],[0,0,0.5,0,0.5],[0,0,0,0,1]] // 2 absorbing states
/*
{ times: [ [ 3 ], [ 3.9999999999999996 ], [ 2.9999999999999996 ] ],
  probs: 
   [ [ 0.75, 0.24999999999999997 ],
     [ 0.49999999999999994, 0.49999999999999994 ],
     [ 0.24999999999999997, 0.7499999999999999 ] ],
  states: [ 1, 5 ] }
*/
			[[0.5,0.25,0.25],[0.5,0,0.5],[0.25,0.25,0.5]]  // no absorbing stats
/*
{ times: [], probs: [], states: [] }
*/
		));
		break;
		
	case "R2.1":  // config methods
		var ran = new RAN({
			trP: [[0.1, 0.9], [0.1, 0.9]]
/*
MRT det= 0.09999999999999998
R> { keys: 
   { index: 'n',
     state: 'u',
     class: 'k',
     x: 'x',
     y: 'y',
     z: 'z',
     t: 't' },
  states: 2,
  syms: { '0': 0, '1': 1 },
  xMap: [ 1, -1 ],
  trProbs: [ [ 0.1, 0.9 ], [ 0.1, 0.9 ] ] }
  */
			//trP: { states: 3, 0: {1: 0.8, 2: 0.1}, 1: {0: 0.1} }
/*
{ 'MRT det': 0 }
Proposed process is not ergodic, thus no unique eq prob exist.  Specify one of the following eq state prs: P^inf -->  [ [ 0.07488979632860272, 0.6664427612064783, 0.25866744246492035 ],
  [ 0.08330534515080978, 0.7413325575350811, 0.1753620973141107 ],
  [ 0, 0, 1 ] ]
R> { keys: 
   { index: 'n',
     state: 'u',
     class: 'k',
     x: 'x',
     y: 'y',
     z: 'z',
     t: 't' },
  states: 3,
  syms: { '0': 0, '1': 1 },
  xMap: [ 0, 1, -1 ],
  trProbs: 
   [ [ 0.09999999999999998, 0.8, 0.1 ],
     [ 0.1, 0.9, 0 ],
     [ 0, 0, 1 ],
     rows: 3,
     columns: 3 ] }
*/
		});
		break;

	case "R2.3":  // config methods
		var ran = new RAN({
			emP: {
				dims: [3,3],
				weights: [1,1]
			},
			trP: { states: 9, 0: {1: 0.8, 2: 0.1}, 1: {0: 0.1}, "0,1": { "1,0": .4} }
		});
/*
trP [ [ 0.09999999999999998, 0.8, 0.1, 0, 0, 0, 0, 0, 0 ],
  [ 0.1, 0.9, 0, 0, 0, 0, 0, 0, 0 ],
  [ 0, 0, 1, 0, 0, 0, 0, 0, 0 ],
  [ 0, 0.4, 0, 0.6, 0, 0, 0, 0, 0 ],
  [ 0, 0, 0, 0, 1, 0, 0, 0, 0 ],
  [ 0, 0, 0, 0, 0, 1, 0, 0, 0 ],
  [ 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
  [ 0, 0, 0, 0, 0, 0, 0, 1, 0 ],
  [ 0, 0, 0, 0, 0, 0, 0, 0, 1 ],
  rows: 9,
  columns: 9 ]
MRT det= 0
Proposed process is not ergodic, thus no unique eq prob exist.  Specify one of the following eq state prs: P^inf -->  [ [ 0.07488979632860272,
    0.6664427612064783,
    0.25866744246492035,
    0,
    0,
    0,
    0,
    0,
    0 ],
  [ 0.08330534515080978,
    0.7413325575350811,
    0.1753620973141107,
    0,
    0,
    0,
    0,
    0,
    0 ],
  [ 0, 0, 1, 0, 0, 0, 0, 0, 0 ],
  [ 0.08545738443447844,
    0.7605083027756313,
    0.15399775120549114,
    0.00003656158440062975,
    0,
    0,
    0,
    0,
    0 ],
  [ 0, 0, 0, 0, 1, 0, 0, 0, 0 ],
  [ 0, 0, 0, 0, 0, 1, 0, 0, 0 ],
  [ 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
  [ 0, 0, 0, 0, 0, 0, 0, 1, 0 ],
  [ 0, 0, 0, 0, 0, 0, 0, 0, 1 ] ]
9 [ [ 0.09999999999999998, 0.8, 0.1, 0, 0, 0, 0, 0, 0 ],
  [ 0.1, 0.9, 0, 0, 0, 0, 0, 0, 0 ],
  [ 0, 0, 1, 0, 0, 0, 0, 0, 0 ],
  [ 0, 0.4, 0, 0.6, 0, 0, 0, 0, 0 ],
  [ 0, 0, 0, 0, 1, 0, 0, 0, 0 ],
  [ 0, 0, 0, 0, 0, 1, 0, 0, 0 ],
  [ 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
  [ 0, 0, 0, 0, 0, 0, 0, 1, 0 ],
  [ 0, 0, 0, 0, 0, 0, 0, 0, 1 ],
  rows: 9,
  columns: 9 ] 1
states [ 0, 1, -1, 2, -2, 3, -3, 4, -4 ]
*/
		break;
		
	case "R2.4":  // config methods
		var ran = new RAN({
			emP: {
				dims: [2,2,2],
				weights: [1,1,1]
			},
			trP: "random"
		});
		break;
		
	case "R3":  // sync pipe with various textbook examples, custom filtering with supervised learning validation
		var ran = new RAN({
			// these have same eqprs [.5, .5] (symmetry -> detailed balance --> eqP[k] = 1/K  eqpr)
			//trP: [[.6, .4],[.4, .6]],
/*
{ stats: 
   { mle_holding_times: [ [Array], [Array], rows: 2, columns: 2 ],
     rel_error: 0.00230801881786998,
     mle_em_probs: null,
     mle_tr_probs: [ [Array], [Array], rows: 2, columns: 2 ],
     tr_counts: [ [Array], [Array], rows: 2, columns: 2 ],
     mean_count: 199.864,
     coherence_time: 5.087382805828901,
     coherence_intervals: 98.28236228402585,
     correlation_0lag: 1,
     mean_intensity: 0.399728,
     degeneracy_param: 2.033569354208375,
     snr: 8.116902388701156 },
  t: 500,
  at: 'end' }
MLEs { holding_Times: '[[0,2.4961510304589924],[2.4806757386729417,0]]',
  viterbiB_emPrs: 'null',
  viterbiA_trPrs: '[[0.601384811290722,0.39861518870927803],[0.40071697978935383,0.5992830202106462]]' }
*/
			//trP: [[0.83177, 0.16822], [0.17152, 0.82848]],
/*
{ stats: 
   { mle_holding_times: [ [Array], [Array], rows: 2, columns: 2 ],
     rel_error: 0.0006768785070861166,
     mle_em_probs: null,
     mle_tr_probs: [ [Array], [Array], rows: 2, columns: 2 ],
     tr_counts: [ [Array], [Array], rows: 2, columns: 2 ],
     mean_count: 84.552,
     coherence_time: 8.834248950434436,
     coherence_intervals: 56.59790694209629,
     correlation_0lag: 1,
     mean_intensity: 0.169104,
     degeneracy_param: 1.493906834514265,
     snr: 5.822665342257233 },
  t: 500,
  at: 'end' }
MLEs { holding_Times: '[[0,5.886354079058032],[5.810975609756097,0]]',
  viterbiB_emPrs: 'null',
  viterbiA_trPrs: '[[0.832333007235839,0.167666992764161],[0.169801047868115,0.830198952131885]]' }
  */
			//trP: [[.5, .5], [.5, .5]],
/*
{ stats: 
   { mle_holding_times: [ [Array], [Array], rows: 2, columns: 2 ],
     rel_error: 0.0006653832163736606,
     mle_em_probs: null,
     mle_tr_probs: [ [Array], [Array], rows: 2, columns: 2 ],
     tr_counts: [ [Array], [Array], rows: 2, columns: 2 ],
     mean_count: 249.94,
     coherence_time: 3.705298163728831,
     coherence_intervals: 134.94190694139024,
     correlation_0lag: 1,
     mean_intensity: 0.49988,
     degeneracy_param: 1.8522044460847678,
     snr: 9.361114481670487 },
  t: 500,
  at: 'end' }
MLEs { holding_Times: '[[0,1.9922586601800394],[1.9936797628135707,0]]',
  viterbiB_emPrs: 'null',
  viterbiA_trPrs: '[[0.5003326916081868,0.4996673083918131],[0.49942810529955756,0.5005718947004424]]' }
  */
			//trP: [[0.1, 0.9], [0.9, 0.1]],

			// textbook exs
			trP: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3
/*  no emP
{ stats: 
   { mle_holding_times: [ [Array], [Array], rows: 2, columns: 2 ],
     rel_error: 0.01865509761388287,
     mle_em_probs: null,
     mle_tr_probs: [ [Array], [Array], rows: 2, columns: 2 ],
     tr_counts: [ [Array], [Array], rows: 2, columns: 2 ],
     mean_count: 91.424,
     coherence_time: 3.4536132327429123,
     coherence_intervals: 144.7759104174188,
     correlation_0lag: 1,
     mean_intensity: 0.182848,
     degeneracy_param: 0.6314862723805761,
     snr: 7.485803061566579 },
  t: 500,
  at: 'end' }
MLEs { holding_Times: '[[0,1.1008115157955753],[9.660927669121593,0]]',
  viterbiB_emPrs: 'null',
  viterbiA_trPrs: '[[0.10186550976138829,0.8981344902386117],[0.10127438873795999,0.89872561126204]]' }
*/
			
			//trP: [[1/2, 1/3, 1/6], [3/4, 0, 1/4], [0,1,0]],  // pg142 ex2  eqpr [.5, .333, .1666]
			//trP: [[1,0,0], [1/4, 1/2, 1/4], [0,0,1]],  // pg143 ex8  no eqprs

			// these have different eqprs
			//trP: [[0.9,0.1],[0.1,0.9]],
			//trP: [[0.1, 0.9], [0.1, 0.9]],  // bernoulli scheme has identical rows
			//trP: [[0.1, 0.9], [0.3, 0.7]],
			//trP: [[0.1, 0.9], [0.4, 0.6]],

			// textbook exs 
			//trP: [[0,1],[1,0]],  // pg433 ex16  regular (all states reachable) absorbing/non on even/odd steps non-regular non-absorbing but ergodic so --> eqpr [.5, .5]
			//trP: [[0.5,0.25,0.25],[0.5,0,0.5],[0.25,0.25,0.5]],  // pg406 ex1  regular (after 2 steps) thus ergodic so eqpr [.4, .2, .4]
			//trP: [[0,1,0,0,0], [0.25,0,0.75,0,0], [0,0.5,0,0.5,0], [0,0,0.75,0,0.25], [0,0,0,1,0]],  // pg433 ex17  non-absorbing non-regular but ergodic so eqpr [.0625, .25, .375]
			//trP: [[1,0,0,0,0],[0.5,0,0.5,0,0],[0,0.5,0,0.5,0],[0,0,0.5,0,0.5],[0,0,0,0,1]],    // 2 absorbing states; non-ergodic so 3 eqpr = [.75 ... .25], [.5 ... .5], [.25 ...  .75]

			//trP: [[1-.2, .1, .1], [0, 1-.1, .1], [.1, .1, 1-.2]],
			//trP: [[1-.2, .1, .1], [0.4, 1-.5, .1], [.1, .1, 1-.2]],
			//trP: [[1-.6, .2, .2,.2], [.1, 1-.3, .1,.1], [.1, .1, 1-.4,.2],[.1,.1,1-.8,.6]],  // non-ergodic
			
			/*
			emP: {
				mu: [ [1], [1.1] ],
				sigma: [ [[1]], [[2]] ]
			}, */
/*
 stats: 
   { mle_holding_times: [ [Array], [Array], rows: 2, columns: 2 ],
     rel_error: 0.015151515151515249,
     mle_em_probs: [ [Object], [Object] ],
     mle_tr_probs: [ [Array], [Array], rows: 2, columns: 2 ],
     tr_counts: [ [Array], [Array], rows: 2, columns: 2 ],
     mean_count: 90.444,
     coherence_time: 3.4312900617373834,
     coherence_intervals: 145.7177886461841,
     correlation_0lag: 1,
     mean_intensity: 0.180888,
     degeneracy_param: 0.6206791966875518,
     snr: 7.470356916777998 },
  t: 500,
  at: 'end' }
MLEs { holding_Times: '[[0,1.0968340824701974],[9.798737174427782,0]]',
  viterbiB_emPrs: '[
  	{"weight":0.5444315148931079,"mu":[1.0139285640388083],"sigma":[[1.020168448272447]],"_gaussian":
  		{"sigma":[[1.020168448272447]],"mu":[1.0139285640388083],"k":1,"_sinv":[[0.9802302763758278]],"_coeff":0.3949791055911031}},
		
	{"weight":0.4555684851068921,"mu":[1.071684787418599],"sigma":[[2.070246038073367]],"_gaussian":
		{"sigma":[[2.070246038073367]],"mu":[1.071684787418599],"k":1,"_sinv":[[0.48303437447011366]],"_coeff":0.2772675754216858}}]',
  viterbiA_trPrs: '[[0.09848484848484848,0.9015151515151515],[0.09996252391565909,0.9000374760843409]]' }
*/

			batch: 50,  // supervised learning every 50 steps
			
			filter: function (str, ev) {  
				switch (ev.at) {
					case "config":
						//Log(ev);
						str.push(ev);
						break;

					case "batch":
						//Log(ev.s,ev.rel_txpr_error);
						Log(ev);
						break;

					case "end":
						Log(ev);
						var
							A = ev.stats.mle_tr_probs,
							B = ev.stats.mle_em_probs,
							H = ev.stats.mle_holding_times;
						
						Log("MLEs", {
							holding_Times: JSON.stringify(H),
							viterbiB_emPrs: JSON.stringify(B),
							viterbiA_trPrs: JSON.stringify(A)
						});
							
						str.push(ev);
						break;
				}
			},

			N: 500,
			steps: 500
		});
		ran.pipe( function (store) {
			//Log(store);
		});
		break;
		
	case "R3.1":  // gen process for R3.2 with async pipe
		var ran = new RAN({

			trP: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3

			batch: 800,  // supervised learning every 50 steps
			
			N: 1000, 
			filter: function (str,ev) {
				switch (ev.at) {
					case "batch":
					case "config":
					case "end":
						Log(JSON.stringify(ev));
				}
			},
			steps: 800  
		});
		ran.pipe(process.stdout);
		/*  
stats : 
{"mle_holding_times":[[0,1.1092640860692102],[9.909690370969315,0]],"rel_error":0.00471644740291818,"mle_em_probs":null,"mle_tr_probs":[[0.10424480266262637,0.8957551973373736],[0.09969867042494697,0.900301329575053]],"tr_counts":[[8394,72128],[71731,647747]],"mean_count":143.859,"coherence_time":5.574154300194505,"coherence_intervals":143.519529047139,"correlation_0lag":1,"mean_intensity":0.17982375,"degeneracy_param":1.0023653293396015,"snr":8.476115384436026}}		
*/
		
		break;		
		
	case "R3.2":  // gen process for R3.3 using async pipe
		var ran = new RAN({

			trP: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3

			//batch: 50,  // supervised learning every 50 steps
			
			N: 10,
			//keys: {state:"u", index: "n"},
			filter: function (str,ev) {
				switch (ev.at) {
					case "jump":
						Log(ev);
						break;
					default:
				}
			},
			steps: 20
		});
		ran.pipe(process.stdout);
		/* copy stdout evs to R3.3 evs */
		break;		
		
	case "R3.3":  // supervised learning with R3.2 evs using asyn pipe
		var 
			getEvents = FLOW.get,	
			evs = [
{ at: 'jump', t: 1, s: 1, index: 3, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 1, s: 1, index: 5, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 1, s: 1, index: 6, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 1, s: 1, index: 7, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 1, s: 1, index: 8, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 1, s: 1, index: 9, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 2, s: 2, index: 2, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 2, s: 2, index: 3, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 3, s: 3, index: 2, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 3, s: 3, index: 6, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 4, s: 4, index: 6, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 5, s: 5, index: 1, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 6, s: 6, index: 1, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 8, s: 8, index: 9, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 9, s: 9, index: 3, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 9, s: 9, index: 8, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 9, s: 9, index: 9, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 10, s: 10, index: 3, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 10, s: 10, index: 8, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 11, s: 11, index: 4, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 13, s: 13, index: 4, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 13, s: 13, index: 8, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 14, s: 14, index: 8, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 15, s: 15, index: 0, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 16, s: 16, index: 0, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 16, s: 16, index: 1, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 16, s: 16, index: 3, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 17, s: 17, index: 3, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 18, s: 18, index: 1, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 18, s: 18, index: 3, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 18, s: 18, index: 6, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 18, s: 18, index: 9, state: 0, hold: 0, obs: null },
{ at: 'jump', t: 19, s: 19, index: 3, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 19, s: 19, index: 6, state: 1, hold: 0, obs: null },
{ at: 'jump', t: 19, s: 19, index: 9, state: 1, hold: 0, obs: null }
				],
			ran = new RAN({

				learn: function (supercb) {
					getEvents(evs, true, (evs) => {
						Log( evs ? ` supervising ${evs.length} events` : " supervised" );
						
						if (evs) // feed supervisor
							supercb(evs);

						else // terminate supervisor
							supercb(null);
					});
				},			

				batch: 50,  // supervised learning every 50 steps

				filter: function (str, ev) {  
					switch (ev.at) {
						case "config":
							Log(ev);
							str.push(ev);
							break;

						case "batch":
							//Log(ev.s,ev.rel_txpr_error);
							Log(ev);
							break;

						case "end":
							//Log(ev);
							str.push(ev);
							break;
					}
				},

				trP: {},  
				//keys: {state:"u", index: "n"},
				K: 2,  // assume 2-state process
				N: 50  // assume 50 members in ensemble
			});
	
		ran.pipe( function (store) {
			//Log(store);
		});
		/*
 1 '00000111110000000000000000000000000000000000000000'
 2 '00010111110000000000000000000000000000000000000000'
 3 '00110101110000000000000000000000000000000000000000'
 4 '00110111110000000000000000000000000000000000000000'
 5 '00110111110000000000000000000000000000000000000000'
 6 '01110111110000000000000000000000000000000000000000'
 8 '01110111100000000000000000000000000000000000000000'
 9 '01100111010000000000000000000000000000000000000000'
10 '01110111110000000000000000000000000000000000000000'
11 '01110111110000000000000000000000000000000000000000'
13 '01111111010000000000000000000000000000000000000000'
14 '01111111110000000000000000000000000000000000000000'
15 '01111111110000000000000000000000000000000000000000'
16 '10101111110000000000000000000000000000000000000000'
17 '10111111110000000000000000000000000000000000000000'
18 '11101101100000000000000000000000000000000000000000'
19 '11111111110000000000000000000000000000000000000000'
*/
		break;	
				
	case "R4.1":  // observation permutations
		Log(perms([],[2,6,4],[]));
		break;
					
	case "R4.2":  // observation permutations
		Log(perms([],[2,6,4],[], function (idx,max) {
			return idx / max;
		}));
		break;
}

// UNCLASSIFIED