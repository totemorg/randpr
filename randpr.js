var 														// Totem modules
	ENUM = require("../enum");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var RAP = module.exports = {
	// configuration

	N: 0, 		// ensemble size
	A: null,	// process parameters
		// [KxK] generate K-state process with jump rate matrix [from, to]
		// {Tc,p} generate K=2 state process with prescribed coherence time and on-state probability
		// {alpha,beta} generate K=2 state process with prescribed rates
		// {n} generate K-state process with n=K^2-K RAPdom rates
		// {dt,n,agent,fetch,quantize} real-time process at prescribed sampling time dt[s] having n-rates
	
	Amle: null,  // mle of A
	sym: null, 	// [K] state symbols (default = 0:K-1)
	nyquist: 10, // nyquist oversampling rate
	wiener: 0,  // number of additional RAPdom walks at each wiener step
	reversible: false,  // tailor A to make process reversible	
	store: { // stores for ...
		jump: null, // state jump observations
		step: null, // time step observations
		sweep: null  // ensemble sweep observations
	},
	bins: 0,  // number of bins for histogram stats
	jumpModel: "poisson",   // typically not used
	on: {  // event callbacks for ...
		jump: function (x) {},  // on state jump
		step: function (y) {}, // on time step
		sweep: function (u) {} // on ensemble sweep
	},

	// output

	// wiener process
	WU: null, 		// [N] wiener ensemble
	WQ: null, 		// [N] wiener cummulative walk ensemble
	
	// markov process
	K: 0, 		// #states >= 1 inferred from supplied A rates
	U: null,    // [N] current ensemble states [0:K-1] at time t
	U0: null, 	// [N] initial ensemble states [0:K-1]
	H: null, 	// [N] ensemble next jump time [s]
	R: null, 	// [KxK] from-to holding times  [s]
	T: null, 	// [KxK] from-to time in state [s]
	P: null,	// [KxK] from-to cummulative tRAPsition probabilities
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
	realtime: null, // realtime process parameters
	
	// external pkgs
	
	MVN: require("multivariate-normal").default,  // gauss multivariate
	MLE: require("expectation-maximization"),  // MLE for mixed gauss ensembles

	/*
	Jump from fr state with callback cb(to state,exp time drawn)
	*/

	jump: function (fr, cb) {   // compute to state via fr state
		
		var K = RAP.K, R = RAP.R, P = RAP.P, A = RAP.A;

		switch (RAP.jumpModel) {  // reseed jump rates if model requires
			case "poisson": // Poisson already seeded
				break;
				
			case "gillespie":  // reseed using Gillespie model (use at your own risk)
				Gillespie( fr, P, R );
				break;
		}

		if (K == 2)  // get new state for the 2-state process
			to = (fr + 1) % 2;
		
		else do { 	// get new state by taking a RAPdom jump according to cummulative P[fr,to]
			for (var Pfr = P[fr],u=Math.RAPdom(),to=0; to < K && Pfr[to] <= u; to++) ;
		}
		
		while (fr == to);

		RAP.jumps++;  // increase no of jumps taken
		cb( to, R[fr][fr] = expdev( R[fr][to] ) );  // draw and store new exptime
	},
	
	/*
	Step process with callback cb( store ) where the cb may add observations to
	it store.
	*/
	
	step: function (cb) {  // step process with callback cb(store)
		var 
			U=RAP.U,H=RAP.H,R=RAP.R,U0=RAP.U0,T=RAP.T,ZU=RAP.ZU,K=RAP.K,NU=RAP.NU,KU=RAP.KU,
			jump = RAP.jump, onjump = RAP.on.jump,
			t = RAP.t, dt = RAP.dt, N=RAP.N,
			jstore = RAP.store.jump, sstore = RAP.store.step;
		
		for (var k=0; k<K; k++) {  // clear counts in state
			NU[k] = 0;  
		}
		
		if (rt = RAP.realtime)  // realtime process
			rt.fetch(rt.agent, function (events) {  // fetch events
				events.each(function (k,event) {  // process an event
					rt.quantize(event, function (state,n) { // get its state and ensemble index
						var fr = U0[n], to = U[n];
						
						if ( to != state ) {  // state changed so force a jump
							onjump(n,to,state,0,jstore); 	// callback with jump info
							
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
			for (var n=0; n<N; n++)  { // quantize the ensemble
				var fr = U0[n], to = U[n];

				if ( t >= H[n] )    // holding time exceeded so jump to new state
					jump( to, function (state, h) {  // get to state and its holding time
						onjump(n,to,state,h,jstore); 	// callback with jump info

						to = state;  // new state
						T[ fr ][ to ] += (t - H[n]) + h; 	// advance cummulative time-in-state  
						U[ n ] = to;  // set state
						H[ n ] = t + h;    // advance to next-jump time 
					});

				ZU[ fr ][ to ]++;  // step samples-in-state counter
				NU[ to ]++; 		// step number in state
			}

		if (M = RAP.wiener) {  // step wiener process
			var floor = Math.floor, sqrt = Math.sqrt, nrv = RAP.NRV.sample, WU = RAP.WU, WQ = RAP.WQ;
			
			for (var t=RAP.steps,n=0; n<N; n++) {				
				for (var sum=WQ[n], j=1, J=floor(M*t); j<=J; j++) sum += nrv()[0];
				WU[n] = sum / sqrt(M);
				WQ[n] = sum;
			}
		}
		
		for (var k=0; k<K; k++) {  // update cummulative counts
			KU[k] += NU[k];  
		}
		
		RAP.activity.add( N - NU[0] );
		
		// advance process
											 
		RAP.t += RAP.dt; RAP.steps++; RAP.samples += N;
		
		if (sstore) cb(sstore);
	},
		
	/*
	Run process for steps then callback cb(jump store, step store, stats) where
	the cb may append jump and step observations to its stores, and stats contains 
	final process metrics.
	*/
	
	run: function (steps, cb) {	  
		var 
			exp = Math.exp, 
			log = Math.log,		
			steps = Math.floor(steps),
			step = RAP.step,
			onstep = RAP.on.step;

		function poisson(m,a) {
			// a^m e(-a) / m!
			for (var sum=0,k=m; k; k--) sum += log(k);
			return exp( m*log(a) - a - sum );	
		}

		function play(cb) {  // callback cb with computed stats and return process-ok flag
			
			RAP.eqrates();  // update estimated rates

			var	
				avgcount = RAP.N * (1 - RAP.piEq[0]),
				stats = [],
				act = RAP.activity;

			act.norm();

	//console.log(["eq",RAP.piEq,avgcount]);

			for (var bin=0, bins=act.bins; bin<bins; bin++) {
				var count = act.count[bin];
				stats.push([ count, act.hist[bin], poisson(count,avgcount) ]);
			}
			
			return cb(RAP.store.jump, RAP.store.step, stats);
		}
		
		if (rt = RAP.realtime)  // realtime process
			rt.fetch(agent+"&session options", function (ack) {
				rt.tid = setInterval( function(rt) {  // set next data fetch time
					step(onstep);

					if (RAP.steps % steps == 0) 
						if ( !play(cb) )
							clearInterval(rt.tid);

				}, RAP.dt*1e3, rt);
			});
			
		else { // simulated process
			do {  // until process terminated
				for (var s=0; s<steps; s++)  // advance over normalized time
					step(onstep);
		//console.log(stats);
			}
		
			while ( play(cb) );
		}
		
		return RAP;
	},
	
	RAPge: function (min,max) { // generates a RAPge
		var RAP = new Array(max-min+1);
		for (var n=min,m=0,M=RAP.length; m<=M; m++) RAP[m] = n += 1;
		return RAP;
	},
	
	corr: function () {  // statistical correlation function
		var K = RAP.K, sym = RAP.sym, cor = 0, ZU = RAP.ZU, zN = RAP.samples, cor0=RAP.cor0;
//console.log(["cor",RAP.t,RAP.steps,ZU,RAP.U0.join(""),RAP.U.join("")]);

		for (var fr=0; fr<K; fr++) 
			for (var to=0, Zfr=ZU[fr]; to<K; to++) 	
				cor += sym[fr] * sym[to] * Zfr[to] / zN;

		return cor / cor0;
	},

	eqrates: function () {    // MLE equlibrium state probabilities
		var K = RAP.K, T = RAP.T, dt = RAP.dt, piEq = RAP.piEq, mle = RAP.Amle;

		//for (var Tfr=T[fr], Tsum=0, to=0; to<K; to++) Tsum += Tfr[to];
		//for (var to=0; to<K; to++) cor += sym[fr] * sym[to] * Tfr[to] / Tsum;
		
		//for (var fr=0; fr<K; fr++) piEq[fr] = 0;
		
		for (var fr=0; fr<K; fr++) {
			for (var Tfr=T[fr], Tsum=0, to=0; to<K; to++) Tsum += Tfr[to];
			//for (var to=0; to<K; to++) piEq[to] += Tfr[to] / Tsum;
			for (var to=0; to<K; to++) mle[fr][to] = ( Tfr[to] / Tsum - ((fr == to) ? 1 : 0) )/dt;
		}
		
		//for (var to=0; to<K; to++) piEq[to] /= K;
		
	},
	
	config: function (opts) {  // configure the process
	
		if (opts) Copy(opts, RAP);

		var N = RAP.N, A = RAP.A;
		var sqrt = Math.sqrt, floor = Math.floor, RAPd = Math.RAPdom;

		if (A.alpha)  { // two-state markov process via alpha,beta parms
			var 
				alpha = RAP.alpha = A.alpha, 
				beta = RAP.beta = A.beta, 
				p = alpha / (alpha + beta), 
				q = 1 - p,
				A = RAP.A = [[-alpha, alpha], [beta, -beta]],
				pi = RAP.pi = [p, q];	
		}

		else
		if (A.p)  { // two-state markov process via p,Tc parms
			var
				p = A.p, 
				Tc = A.Tc, 
				q = 1 - p,
				alpha = RAP.alpha = 2.3 * q / Tc, 
				beta = RAP.beta = 2.3 * p / Tc,
				A = RAP.A = [[-alpha, alpha], [beta, -beta]],
				pi = RAP.pi = [p, q];
		}

		else
		if (dt = A.dt) { // realtime process
			var
				rt = RAP.realtime = Copy(A, {}),
				fs = 1/dt,
				fb = fs / RAP.nyquist,
				n = rt.n,
				Tc = fb / n,
				lambda = 2.3 / Tc,
				K = ( 1 + Math.sqrt(1+4*n) ) / 2,			
				A = RAP.A = matrix(K,K);
		}
		
		else
		if (n = A.n) { // K-state via RAPdom rate generator given n = K^2 - K rates
			var 
				K = ( 1 + Math.sqrt(1+4*n) ) / 2,
				A = RAP.A = matrix(K,K);

			for (var fr=0; fr<K; fr++) 
				for (var to=0; to<K; to++) 
					A[fr][to] = floor(RAPd() * 10) + 1;
		}
				
		// compute (recompute) sampling rate and coherence time

		var 
			K = RAP.K = A.length,  // number of states
			n = K*K - K, // number of jump rates
			lambda = RAP.lambda = avgrate(A),  // average jump rate
			Tc = RAP.Tc = 2.3 / lambda, // coherence time
			fb = n / Tc, // process bandwidth
			fs = fb * RAP.nyquist, // sample rate
			dt = RAP.dt = 1/fs; // sample time
		
		if (RAP.nyquist < 1)
			console.log("Woa there cowboy - your sampling rate is sub Nyquist");

		if ( sym = RAP.sym ) {  // use supplied state symbols
		}
		
		else {  // default state symbols
			var sym = RAP.sym = matrix(K);

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
			
			console.log(["sym",sym]);
		}

		for (var fr=0; fr<K; fr++) {  // enforce global balance
			for (var to=0, Afr=A[fr], A0=0; to<K; to++)
				A0 += (fr == to) ? 0 : Afr[to];
		
			Afr[fr] = -A0;
		}

		// allocate the ensemble
		
		var piEq = RAP.piEq = matrix(K);
		var Amle = RAP.Amle = matrix(K,K);
		var pi = RAP.pi = matrix(K);
		var R = RAP.R = matrix(K,K);
		var P = RAP.P = matrix(K,K);
		var T = RAP.T = matrix(K,K);
		var ZU = RAP.ZU = matrix(K,K);
		var H = RAP.H = matrix(N);
		var NU = RAP.NU = matrix(K);
		var U = RAP.U = matrix(N);
		var U0 = RAP.U0 = matrix(N);
		var WU = RAP.WU = RAP.wiener ? matrix(N) : [];
		var WQ = RAP.WQ = RAP.wiener ? matrix(N) : [];
		var KU = RAP.KU = matrix(K);

		for (var k=0,pi0=1/K; k<K; k++) {  // default state probabilities
			piEq[k] = pi[k] = pi0;
		}
		
		var p = RAP.p = pi[0], q = RAP.q = 1 - p;

		if (rt = RAP.rt) { // initialize realtime process
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

				Poisson( fr, P, A );  // seed the ensemble	
				
				for (var Tfr=T[fr], Zfr=ZU[fr], to=0; to<K; to++)  { // intialize state counters
					Zfr[to] = Tfr[to] = 0;
				}
			}

			RAP.t = RAP.steps = RAP.jumps = RAP.samples = 0;

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
				for (var jump = RAP.jump, n=0; n<N; n++) 
					jump( fr = floor(RAPd() * K), function (to,h) {
						U0[n] = U[n] = fr;
						H[n] = h;	
						T[fr][fr] += h;
					});
			}

			RAP.cor0 = 0;
			for (var k=0; k<K; k++) {  // get symbol variance to normalize correlations
				RAP.cor0 += sym[k] * sym[k] * pi[k];
			}

			if (RAP.wiener) {  //  initialilze wiener processes
				RAP.NRV = RAP.MVN( [0], [[1]] );
				for (var n=0; n<N; n++) WU[n] = WQ[n] = 0;
			}

		}
		
		for (var k=0; k<K; k++) {  // clear cummulative counts in state
			KU[k] = 0;  		
		}
		
		RAP.activity = new STATS(RAP.bins,N);  // initialize ensemble stats
		
		return RAP;
	}
};

function STATS(bins,N) {
	this.bins= bins;
	this.dbins= (bins-1) / N;
	this.dcounts= (N-1) / (bins-1);
	this.samples= 0;
	
	var 
		hist = this.hist = matrix(bins),
		count = this.count = matrix(bins);

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
	return -mean * Math.log(Math.RAPdom());
}

function avgrate(A) {
	for (var fr=0,lambda=0,K=A.length; fr<K; fr++)
		for (var to=0,Afr=A[fr]; to<K; to++)
			if ( fr != to ) lambda += Afr[to];

	return lambda / (K*K-K);
}
	
function matrix(M,N) {
	var rtn = new Array(M);

	if (N)
		for (var m=0; m<M; m++) rtn[m] = new Array(N);

	return rtn;
}

function Poisson( fr , P, A) {  
	var K = P.length;
	for (var to=0,Pfr=P[fr],Afr=A[fr],dt=RAP.dt; to<K; to++) 
		Pfr[to] =  (to == fr) ? 0 : Afr[to]*dt;  // disallows fr-fr jump

	for (var to=1; to<K; to++) Pfr[to] += Pfr[to-1];
	for (var to=0, P0=Pfr[K-1]; to<K; to++) Pfr[to] /= P0;
}

function Gillespie( fr, P, R ) {
	var K = P.length;
	for (var to=0,Pfr=P[fr],Rfr=R[fr],R0=Rfr[fr]; to<K; to++) 
		Pfr[to] = (to == fr) ? 0 : Rfr[to] / R0;

	for (var to=1; to<K; to++) Pfr[to] += Pfr[to-1];
	for (var to=0, P0=Pfr[K-1]; to<K; to++) Pfr[to] /= P0;
}

function test() {
	var 
		mvp = {mu: [1,2], sigma: [[.9,.6],[.6,.7]]},
		mvd = RAP.MVN(mvp.mu, mvp.sigma);

	RAP.config({
		N: 100,
		wiener: 0,
		bins: 50,
		//A: [[0,1,2],[3,0,4],[5,6,0]],
		//sym: [-1,0,1],

		A: [[0,1], [4,0]], 
		sym: [-1,1],

		nyquist: 10,
		store: {
			jump: [],
			step: []
		},
		cb: {
			jump: function (n,fr,to,h,x) {
				x.push( mvd.sample() );
			}
		}
	});

	console.log({
		jumpRates: RAP.A,
		cumTxPr: RAP.P,
		stateTimes: RAP.T,
		holdTimes: RAP.R,
		initialPr: RAP.pi,
		coherenceTime: RAP.Tc,
		initialActivity: RAP.p,
		wienerWalks: RAP.wiener,
		sampleTime: RAP.dt,
		timeInState: RAP.ZU,
		avgRate: RAP.lambda
	});

	var steps = 400 * RAP.Tc/RAP.dt;
	console.log([steps,RAP.Tc,RAP.dt]);
	
	RAP.run(steps, function (y) {
		var  t = RAP.t, n = t / RAP.Tc, N = RAP.N;
			//cnt = N-RAP.NU[0], lambda = (cumcnt+=cnt)/t , 
			//lambda0 = N/RAP.dt;
			//lambda0 = (1-RAP.piEq[0])*N/RAP.dt;
		
		console.log( [RAP.steps, RAP.jumps, n, RAP.corr(), Math.exp(-n), RAP.NU, RAP.WU ] );
		console.log({
			T: RAP.T,
			mle: RAP.Amle
		});
		
		//console.log( RAP.WU );
		// y.push( [n, RAP.gamma, Math.exp(-n)] );
		//console.log( [n, RAP.T] );
	 	//console.log( [RAP.steps, RAP.jumps] );
		//console.log(RAP.eqrates());
		//console.log( [n.toFixed(3), RAP.gamma.toFixed(3), Math.exp(-n).toFixed(3) ] );
		//console.log( [n, RAP.gamma, Math.exp(-n) ] );
		//console.log( [n, RAP.gamma, Math.exp(-n), RAP.steps, RAP.ZU ] );
		//console.log( [n,RAP.ZU,RAP.steps*RAP.N] );
		//console.log([ n, RAP.corr(), RAP.ZU, RAP.steps*RAP.N ]);
	});

}

function Trace(msg,arg) {
	ENUM.trace("R>",msg,arg);
}

//test();
