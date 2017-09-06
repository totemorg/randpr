'use strict';

var 														// Totem modules
	ENUM = require("../enum"),
	STREAM = require("stream");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var LOG = console.log;

class RAN {
	
	constructor(opts) {
		Copy({  // configuration

			N: 0, 		// ensemble size
			A: null,	// process parameters: [dims, units, or range]
				// [KxK] generate K-state process with jump rate MATRIX [from, to]
				// {Tc,p} generate K=2 state process with prescribed coherence time and on-state probability
				// {alpha,beta} generate K=2 state process with prescribed rates
				// {n} generate K-state process with n=K^2-K unqiue random rates
				// {K} generate K-state process with n=K^2-K unqiue random rates
				// {dt,n,agent,fetch,quantize} real-time process at prescribed sampling time dt[s] having n-rates

			Amle: null,  // mle of A: [dims, units, or range]
			sym: null, 	// [K] state symbols (default = 0:K-1)
			nyquist: 10, // nyquist oversampling rate
			wiener: 0,  // number of additional random walks at each wiener step
			reversible: false,  // tailor A to make process reversible	
			bins: 0,  // number of bins for histogram stats
			jumpModel: "poisson",   // typically not used
			intervals: 10, // maximum number of coherence intervals to write streams

			filter: function (str,ev) {  //< streaming filter to modify an event destined for the event stream
				switch ( ev[0] ) {
					case "jump":
						str.push( ev );
						break;
				}
			},

			store: null,   //< provide a [] to make synchronous event stream; null to use async event stream
			
			// outputs

			// wiener process: [dims, units, or range]
			WU: null, 		// [N] wiener ensemble
			WQ: null, 		// [N] wiener cummulative walk ensemble

			// markov process: [dims, units, or range]
			K: 0, 		// #states >= 1 inferred from supplied A rates
			U: null,    // [N] ensemble states [0:K-1] at time t
			U0: null, 	// [N] ensemble states [0:K-1] at time t = 0
			H: null, 	// [N] ensemble next jump time [s]
			R: null, 	// [KxK] from-to holding times  [s]
			T: null, 	// [KxK] from-to time in state [s]
			P: null,	// [KxK] from-to cummulative transition probabilities
			ZU: null, 	// [KxK] from-to samples-in-state probabilities
			pi: null, 	// [K] initial state probabilities (default = 1/K)
			piEq: null, 	// [K] equilibrium  state probabilities
			NU: null, 	// [K] count of ensemble in state [1:N]

			// two-state markov parameters

			alpha: 0,  // on-to-off rate [jumps/s]
			beta: 0,  // off-to-on rate [jumps/s]
			p: 0,  // on state pr 
			q: 0,  // off(not on) state pr 

			// K-state parameters

			lambda: 0,  // average jump rate [jumps/s]
			Tc: 0,  // coherence time >0 [s] 
			dt: 0, // sample time [s]
			t: 0, // step time [s]

			steps: 0, // number of steps 
			jumps: 0, // number of jumps
			samples: 0, // number of elements scanned
			realtime: null // realtime process parameters

		}, this);

		if (opts) this.config(opts);
	}
	
	/*
	Jump from fr state with callback cb(to state,exp time drawn)
	*/
	jump (fr, cb) {   // compute to state via fr state
		
		var 
			K = this.K, R = this.R, P = this.P, A = this.A;

		switch (this.jumpModel) {  // reseed jump rates if model requires
			case "poisson": // Poisson already seeded
				break;
				
			case "gillespie":  // reseed using Gillespie model (use at your own risk)
				Gillespie( this.dt, fr, P, R );
				break;
		}

		if (K == 2)  // get new state for the 2-state process
			to = (fr + 1) % 2;
		
		else do { 	// get new state by taking a random jump according to cummulative P[fr,to]
			for (var Pfr = P[fr],u=Math.random(),to=0; to < K && Pfr[to] <= u; to++) ;
		}
		
		while (fr == to);

		this.jumps++;  // increase no of jumps taken
		cb( to, R[fr][fr] = expdev( R[fr][to] ) );  // draw and store new exptime
	}
	
	/*
	Step process with callback cb( store ) where the cb may add observations to
	it store.
	*/
	step () {  // step process with callback cb(store)
		var 
			ran = this,
			U=this.U,H=this.H,R=this.R,U0=this.U0,T=this.T,ZU=this.ZU,K=this.K,NU=this.NU,KU=this.KU,
			jump = this.jump, onJump = this.onJump,
			t = this.t, dt = this.dt, N=this.N;
		
		for (var k=0; k<K; k++) {  // clear counts in state
			NU[k] = 0;  
		}
		
		/*
		if (rt = this.realtime)  // realtime process
			rt.fetch(rt.agent, function (events) {  // fetch events
				events.each(function (k,event) {  // process an event
					rt.quantize(event, function (state,n) { // get its state and ensemble index
						var fr = U0[n], to = U[n];
						
						if ( to != state ) {  // state changed so force a jump
							onJump(n,to,state,0); 	// callback with jump info
							
							to = state;  // new state
							T[ fr ][ to ] += t - H[n];  // advance cummulative time-in-state  
							U[ n ] = to;	// set state
							H[ n ] = t;	// mark jump time
						}
						
						ZU[ fr ][ to ]++;  // step samples-in-state counter
						NU[ to ]++; 		// step number in state						
					});
				});				
			});  
	
		else // simulated process
		*/
		for (var n=0; n<N; n++)  { // scan the ensemble
			var fr = U0[n], to = U[n];

			if ( t >= H[n] )    // holding time exceeded so jump to new state
				ran.jump( to, function (state, h) {  // get to state and its holding time
					ran.onJump(n,to,state,h); 	// callback with jump info

					to = state;  // new state
					T[ fr ][ to ] += (t - H[n]) + h; 	// advance cummulative time-in-state  
					U[ n ] = to;  		// set state
					H[ n ] = t + h;    // advance to next-jump time 
				});

			ZU[ fr ][ to ]++;  // step samples-in-state counter
			NU[ to ]++; 		// step number-in-state counter
		}

		var M = this.wiener;
		if (M) {  // step wiener process
			var floor = Math.floor, sqrt = Math.sqrt, nrv = this.NRV.sample, WU = this.WU, WQ = this.WQ;
			
			for (var t=this.steps,n=0; n<N; n++) {				
				for (var sum=WQ[n], j=1, J=floor(M*t); j<=J; j++) sum += nrv()[0];
				WU[n] = sum / sqrt(M);
				WQ[n] = sum;
			}
		}
		
		for (var k=0; k<K; k++) {  // update cummulative counts
			KU[k] += NU[k];  
		}
		
		//this.activity.add( N - NU[0] );
		
		// advance process
											 
		this.t += this.dt; this.steps++; this.samples += N;
		
		this.onStep();
	}
		
	/*
	Run process for steps then callback cb(jump store, step store, stats) where
	the cb may append jump and step observations to its stores, and stats contains 
	final process metrics.
	*/
	start (steps, cb) {	  
		var 
			exp = Math.exp, 
			log = Math.log,		
			steps = Math.floor(steps);

		function poisson(m,a) {
			// a^m e(-a) / m!
			for (var sum=0,k=m; k; k--) sum += log(k);
			return exp( m*log(a) - a - sum );	
		}

		/*
		if (rt = this.realtime)  // realtime process
			rt.fetch(agent+"&session options", function (ack) {
				rt.tid = setInterval( function(rt) {  // set next data fetch time
					step();

					if (this.steps % steps == 0) 
						if ( !play(cb) )
							clearInterval(rt.tid);

				}, this.dt*1e3, rt);
			}); 
			
		else { // simulated process
		*/
		
		for (var s=0; s<steps; s++) this.step();  // advance over normalized time

		this.eqrates();  // update estimated rates

		var	
			avgcount = this.N * (1 - this.piEq[0]),
			stats = [],
			act = this.activity;

		act.norm();

//Trace(["eq",this.piEq,avgcount]);

		for (var bin=0, bins=act.bins; bin<bins; bin++) {
			var count = act.count[bin];
			stats.push([ count, act.hist[bin], poisson(count,avgcount) ]);
		}

		this.onStats(stats);
	}
	
	corr () {  // statistical correlation function
		var 
			K = this.K, sym = this.sym, cor = 0, ZU = this.ZU, zN = this.samples, cor0=this.cor0;
//Trace(["cor",this.t,this.steps,ZU,this.U0.join(""),this.U.join("")]);

		for (var fr=0; fr<K; fr++) 
			for (var to=0, Zfr=ZU[fr]; to<K; to++) 	
				cor += sym[fr] * sym[to] * Zfr[to] / zN;

		return cor / cor0;
	}

	eqrates () {    // MLE equlibrium state probabilities
		var 
			K = this.K, T = this.T, dt = this.dt, piEq = this.piEq, mle = this.Amle;

		//for (var Tfr=T[fr], Tsum=0, to=0; to<K; to++) Tsum += Tfr[to];
		//for (var to=0; to<K; to++) cor += sym[fr] * sym[to] * Tfr[to] / Tsum;
		
		//for (var fr=0; fr<K; fr++) piEq[fr] = 0;
		
		for (var fr=0; fr<K; fr++) {
			for (var Tfr=T[fr], Tsum=0, to=0; to<K; to++) Tsum += Tfr[to];
			//for (var to=0; to<K; to++) piEq[to] += Tfr[to] / Tsum;
			for (var to=0; to<K; to++) mle[fr][to] = ( Tfr[to] / Tsum - ((fr == to) ? 1 : 0) )/dt;
		}
		
		//for (var to=0; to<K; to++) piEq[to] /= K;
		
	}

	record (ev) {
		if (this.store) 
			this.filter(this.store, ev);
		
		else
			this.ranStream.push(ev);
	}
	
	onStats (stats) {  //< null to disable
		var n = this.t / this.Tc;
		this.record(["stats",this.t,this.steps, {
			hist: stats,
			number_of_steps: this.steps,
			number_of_state_jumps: this.jumps, 
			coherence_intervals: n, 
			correlation_computed: this.corr(), 
			correlation_theory: Math.exp(-n),
			current_time: this.T,
			jumprate_mle: this.Amle
		}]);
	}

	onJump (idx,to,from,hold) {  // null to disable
		this.record(["jump",this.t,this.steps,idx,to,from,hold]);
	}

	onStep () {		// null to disable
		this.record(["step",this.t, this.steps]);
	}

	pipe(tar, cb) {  // if no cb, stream process events to target stream.  if cb, append transfer events to target list and callback when finished.

		if ( this.store ) {  // in buffering mode
			var steps = this.intervals * this.dt / this.Tc;

			this.start( steps );
			
			if (cb)
				cb( this.store );
			
			else {
				var n=0, N = this.store.length, str = new STEAM.Readable({
					read: function () {
						this.push( this.store[n++] || null );
					}	
				});
				
				str.pipe(tar);
			}
		}

		else 	// in streaming mode
		if (cb) {  // append records to target list
			var tarStream = new STREAM.Writable({
				objectMode: true,
				write: function (ev,en,cb) {
					tar.push(ev);
					cb(null);
				}
			});
				
			this.ranStream.pipe(this.editStream).pipe(tarStream).on("finish", function () {
				cb(tar);
			});
		}
		
		else  // stream records to target stream
			this.ranStream.pipe(this.editStream).pipe(this.charStream).pipe(tar);
	}
		
	onIngest (ev) {

		var
			U=this.U,H=this.H,R=this.R,U0=this.U0,T=this.T,ZU=this.ZU,K=this.K,NU=this.NU,KU=this.KU,		
			n = ev.n,
			t = ev.t,
			fr = U0[n], to = ev.s;

		T[ fr ][ to ] +=  t - H[n]; 	// advance cummulative time-in-state  
		U[ n ] = to;  		// set state
		H[ n ] = t;    		// hold time of last jump

		ZU[ fr ][ to ]++;  // step samples-in-state counter
		NU[ to ]++; 		// step number-in-state counter

	}
	
	config (opts) {  // configure the process
		
		if (opts) Copy(opts, this);
	
		var 
			ran = this,
			N = this.N, A = this.A,
			sqrt = Math.sqrt, floor = Math.floor, rand = Math.random;

		if (A.alpha)  { // two-state markov process via alpha,beta parms
			var 
				alpha = this.alpha = A.alpha, 
				beta = this.beta = A.beta, 
				p = alpha / (alpha + beta), 
				q = 1 - p,
				A = this.A = [[-alpha, alpha], [beta, -beta]],
				pi = this.pi = [p, q];	
		}

		else
		if (A.p)  { // two-state markov process via p,Tc parms
			var
				p = A.p, 
				Tc = A.Tc, 
				q = 1 - p,
				alpha = this.alpha = 2.3 * q / Tc, 
				beta = this.beta = 2.3 * p / Tc,
				A = this.A = [[-alpha, alpha], [beta, -beta]],
				pi = this.pi = [p, q];
		}

		else
		if (dt = A.dt) { // realtime process
			var
				rt = this.realtime = Copy(A, {}),
				fs = 1/dt,
				fb = fs / this.nyquist,
				n = rt.n,
				Tc = fb / n,
				lambda = 2.3 / Tc,
				K = ( 1 + Math.sqrt(1+4*n) ) / 2,			
				A = this.A = MATRIX(K,K);
		}
		
		else
		if (n = A.n) { // K-state via random rate generator given n = K^2 - K rates
			var 
				K = ( 1 + Math.sqrt(1+4*n) ) / 2,
				A = this.A = MATRIX(K,K);

			for (var fr=0; fr<K; fr++) 
				for (var to=0; to<K; to++) 
					A[fr][to] = floor(rand() * 10) + 1;
		}
		else
		if (K = A.K) { // K-state via random rate generator producing n = K^2 - K unique rates
			var 
				n = K*K - K,
				A = this.A = MATRIX(K,K);

			for (var fr=0; fr<K; fr++) 
				for (var to=0; to<K; to++) 
					A[fr][to] = floor(rand() * 10) + 1;
		}
		
		// compute (recompute) sampling rate and coherence time

		var 
			K = this.K = A.length,  // number of states
			n = K*K - K, // number of jump rates
			lambda = this.lambda = avgrate(A),  // average jump rate
			Tc = this.Tc = 2.3 / lambda, // coherence time
			fb = 1 / Tc, // process bandwidth
			fs = fb * this.nyquist, // sample rate
			dt = this.dt = 1/fs; // sample time
		
		//Trace([n,K,lambda,Tc,this.nyquist,dt]);
		
		if (this.nyquist < 1)
			Trace("Woa there cowboy - your sampling rate is sub Nyquist");

		if ( sym = this.sym ) {  // use supplied state symbols
		}
		
		else {  // default state symbols
			var sym = this.sym = MATRIX(K);
			
			if (K % 2) {
				sym[0]=0;
				for (var a=1, k=1; k<K; a++) {
					sym[k++] = a; 
					sym[k++] = -a;
				}
			}
			else			
				for (var a=1, k=0; k<K; a++) {
					sym[k++] = a; 
					sym[k++] = -a;
				}
		}

		for (var fr=0; fr<K; fr++) {  // enforce global balance
			for (var to=0, Afr=A[fr], A0=0; to<K; to++)
				A0 += (fr == to) ? 0 : Afr[to];
		
			Afr[fr] = -A0;
		}

		// allocate the ensemble
		
		var piEq = this.piEq = MATRIX(K);
		var Amle = this.Amle = MATRIX(K,K);
		var pi = this.pi = MATRIX(K);
		var R = this.R = MATRIX(K,K);
		var P = this.P = MATRIX(K,K);
		var T = this.T = MATRIX(K,K);
		var ZU = this.ZU = MATRIX(K,K);
		var H = this.H = MATRIX(N);
		var NU = this.NU = MATRIX(K);
		var U = this.U = MATRIX(N);
		var U0 = this.U0 = MATRIX(N);
		var WU = this.WU = this.wiener ? MATRIX(N) : [];
		var WQ = this.WQ = this.wiener ? MATRIX(N) : [];
		var KU = this.KU = MATRIX(K);

		for (var k=0,pi0=1/K; k<K; k++) {  // default state probabilities
			piEq[k] = pi[k] = pi0;
		}
		
		var p = this.p = pi[0], q = this.q = 1 - p;

		if (rt = this.rt) { // initialize realtime process
			for (var n=0; n<N; n++) {
				H[n] = 0;
				U[n] = U0[n] = -1;
			}
		}
		
		else { // initialize simulated process
			for (var fr=0; fr<K; fr++) {  
				for (var to=0, Afr=A[fr], Rfr=R[fr]; to<K; to++) {  // compute average holding times
					if ( Tc )
						Rfr[to] = (fr == to) ? 0 : 1/Afr[to];

					else
						Rfr[to] = 0;
				}

				Poisson( dt, fr, P, A );  // seed the ensemble	
				
				for (var Tfr=T[fr], Zfr=ZU[fr], to=0; to<K; to++)  { // intialize state counters
					Zfr[to] = Tfr[to] = 0;
				}
			}

			this.t = this.steps = this.jumps = this.samples = 0;

			if (K == 2)  { // initialize two-state process
				for (var n=0, Np=N*p, Ton=R[0][1], Toff=R[1][0]; n<N; n++)  {
					if ( n < Np ) {
						var fr = U0[n] = U[n] = 1;
						var h = H[n] = R[fr][fr] = expdev(Ton);
					}

					else {
						var fr = U0[n] = U[n] = 0;
						var h = H[n] = R[fr][fr] = expdev(Toff);
					}

					T[fr][fr] += h;
				}
			}

			else  { // initialize K-state process
				for (var jump = this.jump, n=0; n<N; n++) 
					jump( fr = floor(rand() * K), function (to,h) {
						U0[n] = U[n] = fr;
						H[n] = h;	
						T[fr][fr] += h;
					});
			}

			this.cor0 = 0;
			for (var k=0; k<K; k++) {  // get symbol variance to normalize correlations
				this.cor0 += sym[k] * sym[k] * pi[k];
			}

			if (this.wiener) {  //  initialilze wiener processes
				this.NRV = RAN.MVN( [0], [[1]] );
				for (var n=0; n<N; n++) WU[n] = WQ[n] = 0;
			}

		}
		
		for (var k=0; k<K; k++) {  // clear cummulative counts in state
			KU[k] = 0;  		
		}
		
		this.activity = new STATS(this.bins,N);  // initialize ensemble stats
		
		//this.ingestStream.on("data", this.onIngest);
		
		// allocate streams
		
		this.ingestStream = new STREAM.Readable({
			objectMode: true
		});

		this.ranStream = new STREAM.Readable({
			objectMode: true,
			read: function () {
				LOG(ran.steps, ran.dt, ran.Tc, ran.intervals);
				
				if ( ran.steps * ran.dt/ran.Tc < ran.intervals ) {			
					//Trace("Buffering 20 steps");
					ran.start(20);
				}
				else
					this.push(null);
			}
		});

		this.editStream = new STREAM.Transform({
			writableObjectMode: true,
			readableObjectMode: true,
			transform: function (ev,en,cb) {
				ran.filter(this,ev);
				cb(null);
			}
		});

		this.charStream = new STREAM.Transform({
			writableObjectMode: true,
			transform: function (ev,en,cb) {
				this.push( JSON.stringify(ev) ); 
				cb(null);
			}
		});

		Trace({  // display config info
			states: this.K,
			ensemble_size: this.N,		
			coherence_interval: this.intervals, 
			coherence_time: this.Tc,
			nyquist_step_time: this.dt,
			jump_rates: this.A,
			cumTxPr: this.P,
			state_times: this.T,
			hold_times: this.R,
			initial_pr: this.pi,
			coherence_time: this.Tc,
			initial_activity: this.p,
			wiener_walks: this.wiener ? "yes" : "no",
			sample_time: this.dt,
			time_in_state: this.ZU,
			avg_rate: this.lambda,
			run_steps: this.intervals * this.dt / this.Tc
		});
	}
	
}

// external pkgs
RAN.MVN = require("multivariate-normal").default; // gauss multivariate
RAN.MLE = require("expectation-maximization");  // MLE for mixed gauss ensembles
			
module.exports = RAN;

function STATS(bins,N) {
	this.bins= bins;
	this.dbins= (bins-1) / N;
	this.dcounts= (N-1) / (bins-1);
	this.samples= 0;
	
	var 
		hist = this.hist = MATRIX(bins),
		count = this.count = MATRIX(bins);

	for (var cnt=1,bin=0; bin<bins; bin++, cnt+=this.dcounts) {
		hist[bin] = 0;
		count[bin] = Math.round(cnt);
	}
}	
	
STATS.prototype.add = function (val) {
	this.hist[ Math.floor( val * this.dbins ) ]++;
	this.samples++;
}

STATS.prototype.norm = function () {
	for (var bin=0, hist=this.hist, bins=this.bins; bin<bins; bin++) hist[bin] /= this.samples * this.dcounts;
}

function expdev(mean) {
	return -mean * Math.log(Math.random());
}

function avgrate(A) {
	for (var fr=0,lambda=0,K=A.length; fr<K; fr++)
		for (var to=0,Afr=A[fr]; to<K; to++)
			if ( fr != to ) lambda += Afr[to];

	return lambda / (K*K-K);
}
	
function MATRIX(M,N) {
	var rtn = new Array(M);

	if (N)
		for (var m=0; m<M; m++) rtn[m] = new Array(N);

	return rtn;
}

function Poisson( dt, fr , P, A) {  
	var K = P.length;
	for (var to=0,Pfr=P[fr],Afr=A[fr]; to<K; to++) 
		Pfr[to] =  (to == fr) ? 0 : Afr[to]*dt;  // disallows fr-fr jump

	for (var to=1; to<K; to++) Pfr[to] += Pfr[to-1];
	for (var to=0, P0=Pfr[K-1]; to<K; to++) Pfr[to] /= P0;
}

function Gillespie( dt, fr, P, R ) {
	var K = P.length;
	for (var to=0,Pfr=P[fr],Rfr=R[fr],R0=Rfr[fr]; to<K; to++) 
		Pfr[to] = (to == fr) ? 0 : Rfr[to] / R0;

	for (var to=1; to<K; to++) Pfr[to] += Pfr[to-1];
	for (var to=0, P0=Pfr[K-1]; to<K; to++) Pfr[to] /= P0;
}

function range (min,max) { // generates a range
	var rtn = new Array(max-min+1);
	for (var n=min,m=0,M=rtn.length; m<=M; m++) rtn[m] = n += 1;
	return rtn;
}	

function test() {
	var 
		ran = new RAN({
			N: 100,
			wiener: 0,
			bins: 50,
			
			//A: [[0,1,2],[3,0,4],[5,6,0]],
			//sym: [-1,0,1],

			A: [[0,1], [4,0]], 
			sym: [-1,1],
			
			store: [],

			intervals: 10,
			nyquist: 4
		});
	
	if (false) 
		ran.pipe(process.stdout);
	
	else
		ran.pipe( [], function (save) {
			LOG(save);
		});
	
}

function Trace(msg,arg) {
	ENUM.trace("R>",msg,arg);
}

//test();
