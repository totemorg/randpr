var 														// Totem modules
	ENUM = require("../enum");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var RAN = module.exports = {
	N: 0, 		// ensemble size
	K: 0, 		// #states
	Ut: null,    // [N] current ensemble states [0:K-1] at time t
	U0: null, 	// [N] initial ensemble states [0:K-1]
	H: null, 	// [N] ensemble next jump times [s]
	R: null, 	// [KxK] from-to holding times  [s]
	T: null, 	// [KxK] from-to jump MLE state times [s]
	A: null,	// [KxK] from-to jump rate matrix [jumps/s]
	P: null,	// [KxK] from-to cummulative transition probabilities
	Z: null, 	// [KxK] from-to state probabilities
	pi: null, 	// [K] initial state probabilities (default = 1/K)
	piEq: null, 	// [K] equilibrium  state probabilities
	E: null, 	// [K] count of ensemble in each state [1:N]
	sym: null, 	// [K] state symbols (default = 0:K-1)

	// two-state markov parameters
	alpha: 0,  // on-to-off rate [jumps/s]
	beta: 0,  // off-to-on rate [jumps/s]
	p: 0,  // on state pr 
	q: 0,  // off(not on) state pr 

	// K-state parameters
	gamma: 1, // ensemble correlation at time t
	lambda: 0,  // average jump rate [jumps/s]
	Tc: 0,  // coherence time [s]
	dt: 0, // sample time [s]
	t: 0, // step time [s]
	nyquist: 10, // nyquist oversampling rate

	jumpModel: "poisson",
	reversible: false,

	steps: 0, // number of steps 
	jumps: 0, // number of jumps
	samples: 0, // number of elements scanned
	
	x: null, // jump observations
	y: null, // step observations
	w: null, // count observations
	
	// external pkgs
	
	MVN: require("multivariate-normal").default,
	MLE: require("expectation-maximization"),
	
	cb: {  // callbacks
		jump: null,  // on jump
		step: null, // after step
		save: null	// after run 
	},

	jump: function (fr, cb) {  // jump from fr state with callback cb(to state,exp time drawn)
		
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
		
		else do { 	// get new state
			for (var Pfr = P[fr],u=Math.random(),to=0; to < K && Pfr[to] <= u; to++) ;
		}
		while (fr == to);

		RAN.jumps++;
		cb( to, R[fr][fr] = expdev( R[fr][to] ) );
	},
	
	step: function (cb) {  // advance the process
		var Ut=RAN.Ut,H=RAN.H,R=RAN.R,U0=RAN.U0,T=RAN.T,Z=RAN.Z,K=RAN.K,E=RAN.E;
		var jump = RAN.jump;
		var t = RAN.t, x = RAN.x, y = RAN.y;
		var jumpcb = x ? RAN.cb.jump : function () {};
		
		for (var k=0; k<K; k++) E[k] = 0;  // initialize counts
		
		for (var n=0, N=RAN.N; n<N; n++)  { // scan the ensemble

			if ( t >= H[n] )    // holding time exceeded so jump to new state
				jump( fr = Ut[n], function (to, h) {  // get new state
					Ut[n] = to;

					jumpcb(n,fr,to,h,x); // callback with jump info
					
					T[ U0[n] ][to] += (t-H[n])+h;  // advance cummulative time-in-state  
					H[n] = t + h;    // advance next-jump time 
				});

			var fr = U0[n], to = Ut[n];
			
			Z[ fr ][ to ]++;  // count the transition
			E[ to ]++; // count number at state
		}

		if (y) cb(y);

		RAN.t += RAN.dt; RAN.steps++; RAN.samples += N;
		RAN.gamma = RAN.corr();
	},
		
	run: function (steps, cb) {	  // run the process for number of steps
		var 
			step = RAN.step,
			stepcb = cb || RAN.cb.step || function () {};
		
		steps = Math.floor(steps);
		while ( steps-- ) {
			step(stepcb);
			if (steps == 10) RAN.eqrates();
		}

		if (save = RAN.cb.save) {
			if (y = RAN.y) save(y,"stepobs");
			if (x = RAN.x) save(x,"jumpobs");
		}

		RAN.eqrates();
		
		return RAN;
	},
	
	range: function (min,max) {
		var ran = new Array(max-min+1);
		for (var n=min,m=0,M=ran.length; m<=M; m++) ran[m] = n += 1;
		return ran;
	},
	
	corr: function () {  // statistical correlation function
		var K = RAN.K, sym = RAN.sym, cor = 0, Z = RAN.Z, zN = RAN.samples, cor0=RAN.cor0;
//console.log(["cor",RAN.t,RAN.steps,Z,RAN.U0.join(""),RAN.Ut.join("")]);

		for (var fr=0; fr<K; fr++) 
			for (var to=0, Zfr=Z[fr]; to<K; to++) 	
				cor += sym[fr] * sym[to] * Zfr[to] / zN;

		return cor / cor0;
	},

	eqrates: function () {    // MLE equlibrium state probabilities
		var K = RAN.K, T = RAN.T, dt = RAN.dt, piEq = RAN.piEq;

		//for (var Tfr=T[fr], Tsum=0, to=0; to<K; to++) Tsum += Tfr[to];
		//for (var to=0; to<K; to++) cor += sym[fr] * sym[to] * Tfr[to] / Tsum;
		
		for (var fr=0; fr<K; fr++) piEq[fr] = 0;
		
		for (var fr=0; fr<K; fr++) {
			for (var Tfr=T[fr], Tsum=0, to=0; to<K; to++) Tsum += Tfr[to];
			for (var to=0; to<K; to++) 	piEq[to] += Tfr[to] / Tsum;
		}
		
		for (var to=0; to<K; to++) piEq[to] /= K;
	},
	
	config: function (opts) {  // configure
	
		if (opts) Copy(opts, RAN);

		var N = RAN.N;

		if (RAN.A) { // K-state markov
		}
		else
		if (RAN.p || RAN.Tc) {  // two-state markov process via p,Tc parms
			var
					p = RAN.p, 
					Tc = RAN.Tc, 
					q = RAN.q = 1 - p,
					alpha = RAN.alpha = 2.3 * q / Tc, 
					beta = RAN.beta = 2.3 * p / Tc;

			RAN.A = [[-alpha, alpha], [beta, -beta]];
			RAN.pi = [p, q];
		}
		else
		if (RAN.alpha || RAN.beta) { // two-state markov process via alpha,beta parms
			var 
				alpha = RAN.alpha, 
				beta = RAN.beta, 
				p = RAN.p = alpha / (alpha + beta), 
				q = RAN.q = 1 - p;

			RAN.A = [[-alpha, alpha], [beta, -beta]];
			RAN.pi = [p, q];	
		}
		else
			console.log("Houston we have a problem");
		
		// compute nyquist sampling rate

		var A = RAN.A, K = RAN.K = A.length;
		var lambda = RAN.lambda = avgrate(A);
		var Tc = RAN.Tc = 2.3 / lambda;
		var dt = RAN.dt = Tc / (K*K-K) / RAN.nyquist;
		
		if (1/dt < 2/Tc)
			console.log("Woa there cowboy - your sampling rate is subNyquist");

		// default state symbols
	
		var sym = RAN.sym;
		if ( !sym ) {
			var sym = matrix(K);
			for (var k=0; k<K; k++) sym[k] = k;
		}

		// enforce global balance
		
		for (var fr=0; fr<K; fr++) {  
			for (var to=0, Afr=A[fr], A0=0; to<K; to++)
				A0 += (fr == to) ? 0 : Afr[to];
		
			Afr[fr] = -A0;
		}

		// allocate the ensemble
		
		var piEq = RAN.piEq = matrix(K);
		var pi = RAN.pi = matrix(K);
		var R = RAN.R = matrix(K,K);
		var P = RAN.P = matrix(K,K);
		var T = RAN.T = matrix(K,K);
		var Z = RAN.Z = matrix(K,K);
		var H = RAN.H = matrix(N);
		var E = RAN.E = matrix(K);
		var Ut = RAN.Ut = matrix(N);
		var U0 = RAN.U0 = matrix(N);

		// default state probabilities
		
		for (var k=0,pi0=1/K; k<K; k++) piEq[k] = pi[k] = pi0;
		
		var p = RAN.p = pi[0], q = RAN.q = 1-p;

		// compute average holding times

		for (var fr=0; fr<K; fr++)   
			for (var to=0, Afr=A[fr], Rfr=R[fr]; to<K; to++) 
				Rfr[to] = (fr == to) ? 0 : 1/Afr[to];
		
		// seed the ensemble

		switch (RAN.jumpModel) {
			case "poisson":
				for (var fr=0; fr<K; fr++) Poisson( fr, P, A );
				break;
				
			case "gillespie":
				for (var fr=0; fr<K; fr++) Gillespie( fr, P, R );
		}
		
		var jump = RAN.jump, floor = Math.floor, rand = Math.random;

		RAN.t = RAN.steps = RAN.jumps = RAN.samples = 0;
		
		for (var fr=0; fr<K; fr++) 
			for (var Tfr=T[fr], Zfr=Z[fr], to=0; to<K; to++) 
				Zfr[to] = Tfr[to] = 0;
		
		if (K == 2)  // two-state 
			for (var n=0, Np=N*p, Ton=R[0][1], Toff=R[1][0]; n<N; n++)  {
				if ( n < Np ) {
					var fr = U0[n] = Ut[n] = 1;
					var h = H[n] = R[fr][fr] = expdev(Ton);
				}
				else {
					var fr = U0[n] = Ut[n] = 0;
					var h = H[n] = R[fr][fr] = expdev(Toff);
				}
				T[fr][fr] += h;
			}

		else  // multi-state 
			for (var n=0; n<N; n++) 
				jump( fr = floor(rand() * K), function (to,h) {
					U0[n] = Ut[n] = fr;
					H[n] = h;	
					T[fr][fr] += h;
				});

		// symbol variance to normalize correlations

		RAN.cor0 = 0;
		for (var k=0; k<K; k++) 
			RAN.cor0 += sym[k] * sym[k] * pi[k];
		
		RAN.gamma = 1;

		return RAN;
	}
};

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
		//A: [[0,1,2],[3,0,4],[5,6,0]],
		//sym: [-1,0,1],

		A: [[0,1],[10,0]], 
		sym: [-1,1],

		nyquist: 10,
		x: [],
		y: [],
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
		initPr: RAN.pi,
		Tc: RAN.Tc,
		p: RAN.p,
		dt: RAN.dt,
		Z: RAN.Z,
		avgRate: RAN.lambda,
		cor: RAN.gamma
	});

	var steps = 4 * RAN.Tc/RAN.dt;
	var cumcnt = 0;
	
	RAN.run(steps, function (y) {
		var  t = RAN.t, n = t / RAN.Tc, N = RAN.N, 
			cnt = N-RAN.E[0], lambda = (cumcnt+=cnt)/t , 
			lambda0 = N/RAN.dt;
			//lambda0 = (1-RAN.piEq[0])*N/RAN.dt;
		
		console.log( [n, RAN.gamma, Math.exp(-n), cnt, lambda / lambda0 ] );
		// y.push( [n, RAN.gamma, Math.exp(-n)] );
		//console.log( [n, RAN.T] );
		//console.log([RAN.steps, RAN.jumps]);
		//console.log(RAN.eqrates());
		//console.log( [n.toFixed(3), RAN.gamma.toFixed(3), Math.exp(-n).toFixed(3) ] );
		//console.log( [n, RAN.gamma, Math.exp(-n) ] );
		//console.log( [n, RAN.gamma, Math.exp(-n), RAN.steps, RAN.Z ] );
		//console.log([n,RAN.Z,RAN.steps*RAN.N]);
		//console.log([ n, RAN.corr(), RAN.Z, RAN.steps*RAN.N ]);
	});

}

//test();
