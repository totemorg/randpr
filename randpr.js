var 														// Totem modules
	ENUM = require("../enum");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var RAN = module.exports = {
	// configuration

	N: 0, 		// ensemble size
	A: null,	// process parameters
		// [KxK] generate K-state process with jump rate matrix [from, to]
		// {Tc,p} generate K=2 state process with prescribed coherence time and on-state probability
		// {alpha,beta} generate K=2 state process with prescribed rates
		// {n} generate K-state process with n=K^2-K random rates
		// {dt,n,agent,fetch,quantize} real-time process at prescribed sampling time dt[s] having n-rates
	
	Amle: null,  // mle of A
	sym: null, 	// [K] state symbols (default = 0:K-1)
	nyquist: 10, // nyquist oversampling rate
	wiener: 0,  // number of additional random walks at each wiener step
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
		
		var K = RAN.K, R = RAN.R, P = RAN.P, A = RAN.A;

		switch (RAN.jumpModel) {  // reseed jump rates if model requires
			case "poisson": // Poisson already seeded
				break;
				
			case "gillespie":  // reseed using Gillespie model (use at your own risk)
				Gillespie( fr, P, R );
				break;
		}

		if (K == 2)  // get new state for the 2-state process
			to = (fr + 1) % 2;
		
		else do { 	// get new state by taking a random jump according to cummulative P[fr,to]
			for (var Pfr = P[fr],u=Math.random(),to=0; to < K && Pfr[to] <= u; to++) ;
		}
		
		while (fr == to);

		RAN.jumps++;  // increase no of jumps taken
		cb( to, R[fr][fr] = expdev( R[fr][to] ) );  // draw and store new exptime
	},
	
	/*
	Step process with callback cb( store ) where the cb may add observations to
	it store.
	*/
	
	step: function (cb) {  // step process with callback cb(store)
		var 
			U=RAN.U,H=RAN.H,R=RAN.R,U0=RAN.U0,T=RAN.T,ZU=RAN.ZU,K=RAN.K,NU=RAN.NU,KU=RAN.KU,
			jump = RAN.jump, onjump = RAN.on.jump,
			t = RAN.t, dt = RAN.dt, N=RAN.N,
			jstore = RAN.store.jump, sstore = RAN.store.step;
		
		for (var k=0; k<K; k++) {  // clear counts in state
			NU[k] = 0;  
		}
		
		if (rt = RAN.realtime)  // realtime process
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

		if (M = RAN.wiener) {  // step wiener process
			var floor = Math.floor, sqrt = Math.sqrt, nrv = RAN.NRV.sample, WU = RAN.WU, WQ = RAN.WQ;
			
			for (var t=RAN.steps,n=0; n<N; n++) {				
				for (var sum=WQ[n], j=1, J=floor(M*t); j<=J; j++) sum += nrv()[0];
				WU[n] = sum / sqrt(M);
				WQ[n] = sum;
			}
		}
		
		for (var k=0; k<K; k++) {  // update cummulative counts
			KU[k] += NU[k];  
		}
		
		RAN.activity.add( N - NU[0] );
		
		// advance process
											 
		RAN.t += RAN.dt; RAN.steps++; RAN.samples += N;
		
		if (sstore) cb(sstore);
	},
		
	/*
	Run process for steps then callback cb(jump store, step store, stats) where
	the cb may append jump and step observations to its stores, and stats contains 
	final process metrics.
	*/
	
	start: function (steps, cb) {	  
		var 
			exp = Math.exp, 
			log = Math.log,		
			steps = Math.floor(steps),
			step = RAN.step,
			onstep = RAN.on.step;

		function poisson(m,a) {
			// a^m e(-a) / m!
			for (var sum=0,k=m; k; k--) sum += log(k);
			return exp( m*log(a) - a - sum );	
		}

		function play(cb) {  // callback cb with computed stats and return process-ok flag
			
			RAN.eqrates();  // update estimated rates

			var	
				avgcount = RAN.N * (1 - RAN.piEq[0]),
				stats = [],
				act = RAN.activity;

			act.norm();

	//console.log(["eq",RAN.piEq,avgcount]);

			for (var bin=0, bins=act.bins; bin<bins; bin++) {
				var count = act.count[bin];
				stats.push([ count, act.hist[bin], poisson(count,avgcount) ]);
			}
			
			return cb(RAN.store.jump, RAN.store.step, stats);
		}
		
		if (rt = RAN.realtime)  // realtime process
			rt.fetch(agent+"&session options", function (ack) {
				rt.tid = setInterval( function(rt) {  // set next data fetch time
					step(onstep);

					if (RAN.steps % steps == 0) 
						if ( !play(cb) )
							clearInterval(rt.tid);

				}, RAN.dt*1e3, rt);
			});
			
		else { // simulated process
			do {  // until process terminated
				for (var s=0; s<steps; s++)  // advance over normalized time
					step(onstep);
		//console.log(stats);
			}
		
			while ( play(cb) );
		}
		
		return RAN;
	},
	
	RAPge: function (min,max) { // generates a RAPge
		var RAN = new Array(max-min+1);
		for (var n=min,m=0,M=RAN.length; m<=M; m++) RAN[m] = n += 1;
		return RAN;
	},
	
	corr: function () {  // statistical correlation function
		var K = RAN.K, sym = RAN.sym, cor = 0, ZU = RAN.ZU, zN = RAN.samples, cor0=RAN.cor0;
//console.log(["cor",RAN.t,RAN.steps,ZU,RAN.U0.join(""),RAN.U.join("")]);

		for (var fr=0; fr<K; fr++) 
			for (var to=0, Zfr=ZU[fr]; to<K; to++) 	
				cor += sym[fr] * sym[to] * Zfr[to] / zN;

		return cor / cor0;
	},

	eqrates: function () {    // MLE equlibrium state probabilities
		var K = RAN.K, T = RAN.T, dt = RAN.dt, piEq = RAN.piEq, mle = RAN.Amle;

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
	
		if (opts) Copy(opts, RAN);

		var N = RAN.N, A = RAN.A;
		var sqrt = Math.sqrt, floor = Math.floor, RAPd = Math.random;

		if (A.alpha)  { // two-state markov process via alpha,beta parms
			var 
				alpha = RAN.alpha = A.alpha, 
				beta = RAN.beta = A.beta, 
				p = alpha / (alpha + beta), 
				q = 1 - p,
				A = RAN.A = [[-alpha, alpha], [beta, -beta]],
				pi = RAN.pi = [p, q];	
		}

		else
		if (A.p)  { // two-state markov process via p,Tc parms
			var
				p = A.p, 
				Tc = A.Tc, 
				q = 1 - p,
				alpha = RAN.alpha = 2.3 * q / Tc, 
				beta = RAN.beta = 2.3 * p / Tc,
				A = RAN.A = [[-alpha, alpha], [beta, -beta]],
				pi = RAN.pi = [p, q];
		}

		else
		if (dt = A.dt) { // realtime process
			var
				rt = RAN.realtime = Copy(A, {}),
				fs = 1/dt,
				fb = fs / RAN.nyquist,
				n = rt.n,
				Tc = fb / n,
				lambda = 2.3 / Tc,
				K = ( 1 + Math.sqrt(1+4*n) ) / 2,			
				A = RAN.A = matrix(K,K);
		}
		
		else
		if (n = A.n) { // K-state via random rate generator given n = K^2 - K rates
			var 
				K = ( 1 + Math.sqrt(1+4*n) ) / 2,
				A = RAN.A = matrix(K,K);

			for (var fr=0; fr<K; fr++) 
				for (var to=0; to<K; to++) 
					A[fr][to] = floor(RAPd() * 10) + 1;
		}
				
		// compute (recompute) sampling rate and coherence time

		var 
			K = RAN.K = A.length,  // number of states
			n = K*K - K, // number of jump rates
			lambda = RAN.lambda = avgrate(A),  // average jump rate
			Tc = RAN.Tc = 2.3 / lambda, // coherence time
			fb = n / Tc, // process bandwidth
			fs = fb * RAN.nyquist, // sample rate
			dt = RAN.dt = 1/fs; // sample time
		
		if (RAN.nyquist < 1)
			console.log("Woa there cowboy - your sampling rate is sub Nyquist");

		if ( sym = RAN.sym ) {  // use supplied state symbols
		}
		
		else {  // default state symbols
			var sym = RAN.sym = matrix(K);

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
		
		var piEq = RAN.piEq = matrix(K);
		var Amle = RAN.Amle = matrix(K,K);
		var pi = RAN.pi = matrix(K);
		var R = RAN.R = matrix(K,K);
		var P = RAN.P = matrix(K,K);
		var T = RAN.T = matrix(K,K);
		var ZU = RAN.ZU = matrix(K,K);
		var H = RAN.H = matrix(N);
		var NU = RAN.NU = matrix(K);
		var U = RAN.U = matrix(N);
		var U0 = RAN.U0 = matrix(N);
		var WU = RAN.WU = RAN.wiener ? matrix(N) : [];
		var WQ = RAN.WQ = RAN.wiener ? matrix(N) : [];
		var KU = RAN.KU = matrix(K);

		for (var k=0,pi0=1/K; k<K; k++) {  // default state probabilities
			piEq[k] = pi[k] = pi0;
		}
		
		var p = RAN.p = pi[0], q = RAN.q = 1 - p;

		if (rt = RAN.rt) { // initialize realtime process
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

			RAN.t = RAN.steps = RAN.jumps = RAN.samples = 0;

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
				for (var jump = RAN.jump, n=0; n<N; n++) 
					jump( fr = floor(RAPd() * K), function (to,h) {
						U0[n] = U[n] = fr;
						H[n] = h;	
						T[fr][fr] += h;
					});
			}

			RAN.cor0 = 0;
			for (var k=0; k<K; k++) {  // get symbol variance to normalize correlations
				RAN.cor0 += sym[k] * sym[k] * pi[k];
			}

			if (RAN.wiener) {  //  initialilze wiener processes
				RAN.NRV = RAN.MVN( [0], [[1]] );
				for (var n=0; n<N; n++) WU[n] = WQ[n] = 0;
			}

		}
		
		for (var k=0; k<K; k++) {  // clear cummulative counts in state
			KU[k] = 0;  		
		}
		
		RAN.activity = new STATS(RAN.bins,N);  // initialize ensemble stats
		
		return RAN;
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
	return -mean * Math.log(Math.random());
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
	for (var to=0,Pfr=P[fr],Afr=A[fr],dt=RAN.dt; to<K; to++) 
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
		mvd = RAN.MVN(mvp.mu, mvp.sigma);

	RAN.config({
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
		jumpRates: RAN.A,
		cumTxPr: RAN.P,
		stateTimes: RAN.T,
		holdTimes: RAN.R,
		initialPr: RAN.pi,
		coherenceTime: RAN.Tc,
		initialActivity: RAN.p,
		wienerWalks: RAN.wiener,
		sampleTime: RAN.dt,
		timeInState: RAN.ZU,
		avgRate: RAN.lambda
	});

	var steps = 400 * RAN.Tc/RAN.dt;
	console.log([steps,RAN.Tc,RAN.dt]);
	
	RAN.start(steps, function (y) {
		var  t = RAN.t, n = t / RAN.Tc, N = RAN.N;
			//cnt = N-RAN.NU[0], lambda = (cumcnt+=cnt)/t , 
			//lambda0 = N/RAN.dt;
			//lambda0 = (1-RAN.piEq[0])*N/RAN.dt;
		
		console.log( [RAN.steps, RAN.jumps, n, RAN.corr(), Math.exp(-n), RAN.NU, RAN.WU ] );
		console.log({
			T: RAN.T,
			mle: RAN.Amle
		});
		
		//console.log( RAN.WU );
		// y.push( [n, RAN.gamma, Math.exp(-n)] );
		//console.log( [n, RAN.T] );
	 	//console.log( [RAN.steps, RAN.jumps] );
		//console.log(RAN.eqrates());
		//console.log( [n.toFixed(3), RAN.gamma.toFixed(3), Math.exp(-n).toFixed(3) ] );
		//console.log( [n, RAN.gamma, Math.exp(-n) ] );
		//console.log( [n, RAN.gamma, Math.exp(-n), RAN.steps, RAN.ZU ] );
		//console.log( [n,RAN.ZU,RAN.steps*RAN.N] );
		//console.log([ n, RAN.corr(), RAN.ZU, RAN.steps*RAN.N ]);
	});

}

function Trace(msg,arg) {
	ENUM.trace("R>",msg,arg);
}

test();
