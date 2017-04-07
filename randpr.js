var 														// Totem modules
	ENUM = require("../enum");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var RAN = module.exports = {
	N: 0, 		// ensemble size
	K: 0, 		// #states
	E: null,    // [N] current ensemble states [0:K-1]
	I: null, 	// [N] initial ensemble states [0:K-1]
	H: null, 	// [N] ensemble next jump times [s]
	R: null, 	// [KxK] from-to holding times  [s]
	T: null, 	// [KxK] from-to jump state times [s]
	A: null,	// [KxK] from-to jump rate matrix [jumps/s]
	W: null, 	// [KxK] from-to state transition probabilities
	P: null,	// [KxK] from-to cummulative transition pr matrix
	pi: null, 	// [K] equilibrium probs (default = 1/K)
	sym: null, 	// [K] state symbols (default = 0:K-1)

	// two-state markov parameters
	alpha: 0,  // on-to-off rate [jumps/s]
	beta: 0,  // off-to-on rate [jumps/s]
	p: 0,  // on state pr 
	q: 0,  // off(not on) state pr 

	// K-state parameters
	lambda: 0,  // average jump rate [jumps/s]
	Tc: 0,  // coherence time [s]
	dt: 0, // sample time [s]
	t: 0, // step time [s]
	nyquist: 10, // nyquist oversampling rate

	jumpModel: "poisson",
	steps: 0, // step counter	
	jumps: 0, // number of jumps
	reversible: false,
	x: null, // jump observations
	y: null, // step observations

	// external pkgs
	
	MVN: require("multivariate-normal").default,
	MLE: require("expectation-maximization"),
	
	cb: {  // callbacks
		jump: null,  // on jump
		step: null, 	// on step
		save: null		// on run end
	},

	jump: function (fr, cb) {  // jump from fr state with callback cb(to state,exp time drawn)
		
		var K = RAN.K, R = RAN.R, P = RAN.P, A = RAN.A;

		switch (RAN.jumpModel) {
			case "poisson": // Poisson jump model already initialized on config
				break;
				
			case "gillispie":  // Gillispie jump model (use at your own risk)
				Gillispie( fr, P, R );
				break;
		}

		if (K == 2)  // get new state
			to = (fr + 1) % 2;
		
		else do { 	// get new state
			for (var Pfr = P[fr],u=Math.random(),to=0; to < K && Pfr[to] <= u; to++) ;
		}
		while (fr == to);

		RAN.jumps++;
		cb( to, R[fr][fr] = expdev( R[fr][to] ) );
		return to;
	},
	
	step: function (cb) {  // advance the process
		var E=RAN.E,H=RAN.H,R=RAN.R,I=RAN.I,T=RAN.T,Z=RAN.Z;
		var jump=RAN.jump;
		var t = RAN.t; // += RAN.dt, steps = RAN.steps++;
		var x = RAN.x, y = RAN.y;
		var jumpcb = x ? RAN.cb.jump : function () {};
		
		for (var n=0, N=RAN.N; n<N; n++)  { // scan the ensemble

			Z[ I[n] ][ E[n] ]++;

			if ( t >= H[n] )    // holding time exceeded so jump to new state
				var to = jump( fr = E[n], function (to, h) {  // get new state
					E[n] = to;

					jumpcb(n,fr,to,h,x); // callback with jump info
					
					//T[ I[n] ][to] += h;  // advance cummulative time-in-state  
					//H[n] = h;    // advance next-jump time 
					//T[ fr ][to] = (t-H[n])+h;  // advance cummulative time-in-state  
					T[ I[n] ][to] += (t-H[n])+h;  // advance cummulative time-in-state  
					H[n] = t + h;    // advance next-jump time 
					//Z[ I[n] ][to]++;
				});
		}							

		if (y) cb(y);

		RAN.t += RAN.dt; RAN.steps++;
		RAN.gamma = RAN.corr();

		//RAN.t += RAN.dt;  RAN.steps++;
	},
		
	run: function (steps, cb) {	  // run the process for number of steps
		var 
			step = RAN.step,
			stepcb = cb || RAN.cb.step || function () {};
		
		steps = Math.floor(steps);
		while ( steps-- ) step(stepcb);

		if (save = RAN.cb.save) {
			if (y = RAN.y) save(y,"stepobs");
			if (x = RAN.x) save(x,"jumpobs");
		}

		return RAN;
	},
	
	corr: function () {  // statistical correlation function
		var K = RAN.K, T = RAN.T, sym = RAN.sym, cor = 0, corr0 = RAN.corr0 , Z = RAN.Z;

		var zN = RAN.steps* RAN.N;
//console.log(["cor",RAN.t,RAN.steps,Z,RAN.I.join(""),RAN.E.join("")]);

		for (var fr=0; fr<K; fr++) {
			//for (var Tfr=T[fr], Tsum=0, to=0; to<K; to++) Tsum += Tfr[to];
			//for (var to=0; to<K; to++) 	cor += sym[fr] * sym[to] * Tfr[to] / Tsum;
			for (var to=0, Zfr=Z[fr]; to<K; to++) 	cor += sym[fr] * sym[to] * Zfr[to] / zN;
		}

		return cor ; // / corr0 ;
	},

	txprs: function () {    // compute state transition probs
		var K = RAN.K, T = RAN.T, W = RAN.W;

		for (var fr=0; fr<K; fr++) {
			for (var Tfr=T[fr], Wfr=W[fr], Tsum=0, to=0; to<K; to++) Tsum += Tfr[to];
			for (var to=0; to<K; to++) 	Wfr[to] = Tfr[to] / Tsum;
		}

		return RAN.W;		
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

		// default equilib probabilities
		
		var pi = RAN.pi;		
		if ( !pi ) {
			var pi = RAN.pi = matrix(K);
			for (var fr=0,pi0=1/K; fr<K; fr++) pi[fr] = pi0;
		}
		
		var p = RAN.p = pi[0], q = RAN.q = 1-p;

		// enforce global balance
		
		for (var fr=0; fr<K; fr++) {  
			for (var to=0, Afr=A[fr], A0=0; to<K; to++)
				A0 += (fr == to) ? 0 : Afr[to];
		
			Afr[fr] = -A0;
		}

		// allocate the ensemble
		
		var R = RAN.R = matrix(K,K);		
		var P = RAN.P = matrix(K,K);
		var T = RAN.T = matrix(K,K);
		var W = RAN.W = matrix(K,K);
		var E = RAN.E = matrix(N);
		var H = RAN.H = matrix(N);
		var I = RAN.I = matrix(N);
		var Z = RAN.Z = matrix(K,K);

		// compute average holding times

		for (var fr=0; fr<K; fr++)   
			for (var to=0, Afr=A[fr], Rfr=R[fr]; to<K; to++) 
				Rfr[to] = (fr == to) ? 0 : 1/Afr[to];
		
		// seed the ensemble

		for (var fr=0; fr<K; fr++) Poisson( fr, P, A );

		var jump = RAN.jump, floor = Math.floor, rand = Math.random;

		RAN.t = RAN.steps = RAN.jumps = 0;
		
		for (var fr=0; fr<K; fr++) for (var Tfr=T[fr], Zfr=Z[fr], to=0; to<K; to++) Zfr[to] = Tfr[to] = 0;
		
		if (K == 2)  // two-state 
			for (var n=0, Np=N*p, Ton=R[0][1], Toff=R[1][0]; n<N; n++)  {
				if ( n < Np ) {
					var fr = I[n] = to = E[n] = 1;
					var h = H[n] = R[fr][fr] = expdev(Ton);
				}
				else {
					var fr = I[n] = to = E[n] = 0;
					var h = H[n] = R[fr][fr] = expdev(Toff);
				}
				T[fr][fr] += h;
				//Z[fr][fr] ++;
			}

		else  // multi-state 
			for (var n=0; n<N; n++) 
				jump( fr = floor(rand() * K), function (to,h) {
					I[n] = fr;
					E[n] = to;
					H[n] = h;	
					//T[fr][fr] += h;
				});

		// initial symbol-value normalized correlation

		RAN.corr0 = 0;
		for (var k=0; k<K; k++) RAN.corr0 += sym[k] * sym[k];
			
		//RAN.steps = 0; RAN.t = 0;
		RAN.gamma = 1;
		//RAN.run(1); RAN.t  = RAN.steps = RAN.jumps = 0;

		RAN.txprs();

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

function Gillispie( fr, P, R ) {
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
		N: 10,
		//A: [[0,1,2],[3,0,4],[5,6,0]],
		//sym: [-1,0,1],

		A: [[0,1],[1,0]], 
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
		eqPr: RAN.pi,
		Tc: RAN.Tc,
		p: RAN.p,
		dt: RAN.dt,
		Z: RAN.Z,
		cor: RAN.gamma
	});

	var steps = 4 * RAN.Tc/RAN.dt;

	RAN.run(steps, function (y) {
		var n = RAN.t / RAN.Tc;
		// y.push( [n, RAN.gamma, Math.exp(-n)] );
		//console.log( [n, RAN.T] );
		//console.log([RAN.steps, RAN.jumps]);
		//console.log(RAN.txprs());
		//console.log( [n.toFixed(3), RAN.gamma.toFixed(3), Math.exp(-n).toFixed(3) ] );
		//console.log( [n, RAN.gamma, Math.exp(-n) ] );
		console.log( [n, RAN.gamma, Math.exp(-n), RAN.steps, RAN.Z ] );
		//console.log([n,RAN.Z,RAN.steps*RAN.N]);
		//console.log([ n, RAN.corr(), RAN.Z, RAN.steps*RAN.N ]);
	});

}

//stest();


