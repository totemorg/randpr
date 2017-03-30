var 														// Totem modules
	ENUM = require("../enum");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var RP = module.exports = {
	E: [],    // ensemble 
	H: [], 	// ensemble holding times [s]
	N: 50, // ensemble size
	K: 2, // #states

	// two-state markov parameters
	alpha: 0,  // on-to-off rate [jumps/s]
	beta: 0,  // off-to-on rate [jumps/s]
	lambda: 0,  // jump rate [jumps/s]
	p: 0,  // on state pr (activity)
	q: 0,  // off state pr (inactivity)
	Tc: 0,  // coherence time [s]
	dT: 0, // sample time [s]
	
	// multi-state gauss mixing parameters
	mu: null,  // Mx1 mean lists [m]
	sigma: null,  // Mx1 sigma list [m]
	A: null,  // KxK from-to jump rate matrix [jumps/s]
	P: null,  // KxK from-to cummulative transition pr matrix
	
	gill: function (fr, A, P) {

		if (K = A.length == 2)
			return (fr + 1) % 2;

		for (var u=Math.rand(),k=0; k<K; k++) if ( P[k] > u ) return k;
		return fr;
	},

	config: function (opts, cb) {
	
		if (opts) Copy(opts, RP);

		var N = RP.N;
		var E = RP.E = new Array(N);
		var H = RP.H = new Array(N);
		var P = RP.P = new Array(K);
		
		if (RP.mu || RP.sigma) { // multi-state guass mixing process
			var A = RP.A = [];
		}
		else
		if (p = RP.p || Tc = RP.Tc) {  // two-state markov process via p,Tc parms
			var q = RP.q = 1 - p;
			var K = RP.K = 2;
			var alpha = RP.alpha = 2.3 * q / Tc;
			var beta = RP.beta = 2.3 * p / Tc;
			var lambda = RP.lambda = (alpha + beta) / 2;
			var A = RP.A = [[1-alpha, alpha], [beta, 1-beta]];
		}
		else
		if (alpha = RP.alpha || beta = TP.beta) { // two-state markov process via alpha,beta parms
			var K = RP.K = 2;			
			var p = RP.p = alpha / (alpha + beta);
			var q = RP.q = beta / (alpha + beta);
			var Tc = RP.Tc = 2.3 / lambda;
			var lambda = RP.lambda = (alpha + beta) / 2;
			var A = RP.A = [[1-alpha, alpha], [beta, 1-beta]];
		}
		else
			console.log("Huston we have a problem");
		
		if (1/RP.dT < 2/RP.Tc)
			console.log("Woo there cowboy - your frame rate is subNyquist");
		
		// define cummulative from-to pr transition matrix using Gillispie model

		var P = RP.P = new Array(K);
		for (var k=0; k<K; k++) P[k] = new Array(K);

		for (var fr=0; fr<K; fr++) 
			for (var B0=0,to=0,Pfr=P[fr],Afr=A[fr]; to<K; B0+=Pfr[to], to++)
				Pfr[to] = B0 + (to == fr) ? 0 : Afr[to] / Afr[fr];
	
		for (var to=0; to<K; to++) Pfr[to] /= B0;

		// initialize the ensemble

		if (mu || sigma) // multi-state guass mix
			for (var gill=RP.gill, n=0; n<N; n++) {
				var 
					fr = Math.floor(Math.rand() * K),
					to = gill( fr, A[fr], B );
				
				E[n] = fr;
				H[n] = expdev(1/A[fr][to]);
			}
					
		else  // two-state 
			for (var n=0, Np=N*p, Ton=1/alpha, Toff=1/beta; n<N; n++) 
				if ( n <= Np ) {
					E[n] = 1;
					H[n] = expdev(Ton);
				}
				else {
					E[n] = 0;
					H[n] = expdev(Toff);
				}
		
		console.log(RP);
	}
};

function expdev(mean) {
	return -mean * Math.log(Math.random());
}

RP.config({
	N: 10,
	p: 0.1,
	Tc: 1
});
