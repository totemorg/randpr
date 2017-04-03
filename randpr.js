var 														// Totem modules
	ENUM = require("../enum"),
	MVN = require("multivariate-normal").default,
	EM = require("expectation-maximization");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var RP = module.exports = {
	E: null,    // N ensemble 
	H: null, 	// N ensemble exit times [s]
	R: null, 	// KxK from-to holding times  [s]
	N: 0, 		// ensemble size
	K: 0, 		// #states
	J: null, 	// KxK from-to jump counts
	A: null,	// KxK from-to jump rate matrix [jumps/s]
	P: null,	// KxK from-to cummulative transition pr matrix
	pi: null, 	// K equilibrium probs

	// two-state markov parameters
	alpha: 0,  // on-to-off rate [jumps/s]
	beta: 0,  // off-to-on rate [jumps/s]
	lambda: 0,  // jump rate [jumps/s]
	p: 0,  // on state pr (activity)
	q: 0,  // off state pr (inactivity)
	Tc: 1000,  // coherence time [s]
	dt: 1, // sample time [s]
	t: 0, // step time [s]
	jumpModel: "poisson",
	nyquist: 10, // nyquist oversampling rate
	steps: 0, // step counter	
	jumps: 0, // number of jumps
	mix: null,  // multi-state mixing parameters
	x: null, // mix observations
	reversible: false,
	
	jump: function (fr) {  
		
		var K = RP.K, R = RP.R, P = RP.P, A = RP.A;

		function Poisson( fr ) {
			for (var to=0,Pfr=P[fr],Afr=A[fr],dt=RP.dt; to<K; to++) 
				Pfr[to] =  (to == fr) ? 0 : Afr[to]*dt;  // disallows fr-fr jump
		
			//console.log(["pb",fr,Pfr]);
			for (var to=1; to<K; to++) Pfr[to] += Pfr[to-1];
			for (var to=0, P0=Pfr[K-1]; to<K; to++) Pfr[to] /= P0;
			//console.log(["pa",fr,Pfr]);
		}
		
		function Gillispie( fr ) {
			for (var to=0,Pfr=P[fr],Rfr=R[fr],R0=Rfr[fr]; to<K; to++) 
				Pfr[to] = (to == fr) ? 0 : Rfr[to] / R0;

			//console.log(["gb",fr,Pfr]);
			for (var to=1; to<K; to++) Pfr[to] += Pfr[to-1];
			for (var to=0, P0=Pfr[K-1]; to<K; to++) Pfr[to] /= P0;
			//console.log(["ga",fr,Pfr]);
		}
		
		if (K == 2) return (fr + 1) % 2;

		switch (RP.jumpModel) {
			case "off":
				break;
				
			case "poisson": // Poisson jump model
				Poisson( fr );
				break;
				
			case "gillispie":  // Gillispie jump model
				Gillispie( fr );
				break;
		}

		for (var Pfr = P[fr],u=Math.random(),to=0; to < K && Pfr[to] <= u; to++) ;

		R[fr][fr] = expdev(R[fr][to]);
		RP.jumps++;
		RP.J[fr][to]++;
		
		return to;
	},
	
	step: function (cb) {
		var E=RP.E,H=RP.H,R=RP.R,jump=RP.jump;
		var t = RP.t += RP.dt, steps = RP.steps++;
		var x = RP.x;
		
		for (var n=0,N=RP.N; n<N; n++) {
			if ( t > H[n] ) {
				var fr = E[n], to = E[n] = jump( fr ), h = R[fr][to];
				H[n] += expdev( h );
				if (cb) cb(n,fr,to,h,x);
				//console.log([t,n,fr,to,H[n]]);
			}
			//else
				//console.log([t,n,E[n],E[n],H[n]]);
		}
		//console.log([RP.steps,E]);
	},
		
	run: function (steps, cb) {	
		var step = RP.step;
		
		while (steps--) step(cb);
		return RP;
	},
	
	reset: function () {
		RP.t = 0;
		return RP;
	},
	
	meanRate: function (A) {
		for (var fr=0,lambda=0,K=RP.K; fr<K; fr++)
			for (var to=0, Afr=A[fr]; to<K; to++)
				if ( fr != to ) lambda += Afr[to];
		
		return lambda/(K*K-K);
	},
	
	config: function (opts, cb) {
	
		function matrix(M,N) {
			var rtn = new Array(M);
			if (N)
				for (var m=0; m<M; m++) rtn[m] = new Array(N);
			return rtn;
		}

		if (opts) Copy(opts, RP);

		var N = RP.N;
		var E = RP.E = matrix(N);
		var H = RP.H = matrix(N);

		if (RP.A) { // K-state markov
			var A = RP.A, K = RP.K = A.length;
			var lambda = RP.lambda = RP.meanRate(A);
			var Tc = RP.Tc = 2.3/lambda;
		}
		else
		if (RP.p || RP.Tc) {  // two-state markov process via p,Tc parms
			var p = RP.p, Tc = RP.Tc, q = RP.q = 1 - p, K = RP.K = 2, dt = RP.dt;
			var alpha = RP.alpha = 2.3 * q / Tc, beta = RP.beta = 2.3 * p / Tc, lambda = RP.lambda = (alpha + beta) / 2;
			var A = RP.A = [[-alpha, alpha], [beta, -beta]];
		}
		else
		if (RP.alpha || RP.beta) { // two-state markov process via alpha,beta parms
			var alpha = RP.alpha, beta = RP.beta, K = RP.K = 2, dt = RP.dt;
			var p = RP.p = alpha / (alpha + beta), q = RP.q = 1 - p;
			var Tc = RP.Tc = 2.3 / lambda, lambda = RP.lambda = (alpha + beta) / 2;
			var A = RP.A = [[-alpha, alpha], [beta, -beta]];
		}
		else
			console.log("Huston we have a problem");
		
		var dt = RP.dt = Tc/RP.nyquist;
		
		if (1/dt < 2/Tc)
			console.log("Woa there cowboy - your sampling rate is subNyquist");

		var P = RP.P = matrix(K,K);
		var R = RP.R = matrix(K,K);
		var J = RP.J = matrix(K,K);
		
		// default equilib probabilities
		
		var pi = RP.pi;
		
		if ( !pi ) {
			var pi = RP.pi = matrix(K);
			for (var fr=0,pi0=1/K; fr<K; fr++) pi[fr] = pi0;
		}
		
		console.log(["pi",pi]);
		
		// enforce global balance
		
		for (var fr=0; fr<K; fr++) {  
			for (var to=0, Afr=A[fr], A0=0; to<K; to++)
				A0 += (fr == to) ? 0 : Afr[to];
		
			Afr[fr] = -A0;
		}

		console.log(["Ainit",A]);
		
		// compute holding times
		
		for (var fr=0; fr<K; fr++)   
			for (var to=0, Afr=A[fr], Rfr=R[fr]; to<K; to++) 
				Rfr[to] = (fr == to) ? 0 : 1/Afr[to];
		
		console.log(["Rinit", Tc, R]);

		// initialize poisson jump model
		var jump=RP.jump;
		
		for (var fr=0; fr<K; fr++) jump(fr);
		RP.jumpModel = "off";
		
		console.log( ["Pinit",dt,P] );
		
		// seed the ensemble
		var floor = Math.floor, rand=Math.random;
		
		if (K == 2)  // two-state 
			for (var n=0, Np=N*p, Ton=R[0][1], Toff=R[1][0]; n<N; n++)  {
				if ( n <= Np ) {
					E[n] = 1;
					H[n] = expdev(Ton);
				}
				else {
					E[n] = 0;
					H[n] = expdev(Toff);
				}
			}

		else  // multi-state guass mix
			for (var n=0; n<N; n++) {
				E[n] = jump( fr = floor(rand() * K) );
				H[n] = R[fr][fr];				
			}

		// initialize jump counters
		
		for (var fr=0; fr<K; fr++) for (var Jfr=J[fr],to=0; to<K; to++) 
			Jfr[to] = (fr == to) ? floor( pi[fr] * K ) : 0;
		
		console.log(["Jinit",J]);
		
		RP.jumps = K;
		
		return RP;
	}
};

function expdev(mean) {
	return -mean * Math.log(Math.random());
}

var 
	mvp = {mu: [1,2], sigma: [[.9,.6],[.6,.7]]},
	mvd = MVN(mvp.mu, mvp.sigma);

RP.config({
	N: 5,
	A: [[0,1,2],[3,0,4],[5,6,0]],
	nyquist: 20,
	mix: mvp,
	x: []
	//p: 0,
	//Tc: 1
}).run(8, function (n,fr,to,h,x) {
	x.push( mvd.sample() );
	console.log(["jump",RP.jumps,RP.steps,n,fr,to,h,RP.J]);	
});

