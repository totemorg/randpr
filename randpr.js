var 														// Totem modules
	ENUM = require("../enum");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var RP = module.exports = {
	E: [],    // ensemble 
	H: [], 	// holding times [s]
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
	mu: null,  // means   [m]
	sigma: null,  // sigmas  [m]
	A: null,  .. // jump rate matrix [from, ... ][to, ...] in [jumps/s]
	B: null,
	
	gill: function (fr, A, B) {
		for (var B0=0,k=0,K=A.length; k<K; B0+=B[k], k++)
			B[k] = B0 + (k == fr) ? 0 : A[k] / A[fr];	
	
		for (var k=0; k<K; k++) B[k] /= B0;
	
		for (var u=Math.rand(),k=0; k<K; k++) if ( B[k] > u ) return k;
		return fr;
	},

	config: function (opts, cb) {
	
		if (opts) Copy(opts, RP);

		var N = RP.N;
		var E = RP.E = new Array(N);
		var H = RP.H = new Array(N);
		var B = RP.B = new Array(K);
		
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