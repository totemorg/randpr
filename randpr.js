var 														// Totem modules
	ENUM = require("../enum");
	
var															// shortcuts
	Copy = ENUM.copy,
	Each = ENUM.each;

var RP = module.exports = {
	E: [],    // ensemble
	H: [], 	// holding times
	N: 50,
	alpha: 0,  // on-to-off jump rate
	beta: 0,  // off-to-on jump rate
	lambda: 0,  
	p: 0,  // on state pr
	q: 0,  // off state pr
	Tc: 0,  // coherence time
	dT: 0, // sample time
	config: function (opts) {
		if (opts) Copy(opts, RP);

		RP.E = new Array(RP.N);
		RP.H = new Array(RP.N);
		
		if (RP.p || RP.Tc) {
			RP.q = 1 - RP.p;
			RP.alpha = 2.3 * RP.q / RP.Tc;
			RP.beta   = 2.3 * RP.p / RP.Tc;
		}
		else {
			RP.p = RP.alpha / (RP.alpha + RP.beta);
			RP.q = RP.beta / (RP.alpha + RP.beta);
			RP.Tc = 2.3 / RP.lambda;
		}
		
		for (var n=0, N=RP.N, Np=RP.N*RP.p, Ton=1/RP.alpha, Toff=1/RP.beta; n<N; n++) 
			if ( n <= Np ) {
				RP.E[n] =1;
				RP.H[n] = expdev(Ton);
			}
			else {
				RP.E[n] = 0;
				RP.H[n] = expdev(Toff);
			}
					
		RP.lambda = (RP.alpha + RP.beta) / 2;
		
		if (1/RP.dT < 2/RP.Tc)
			console.log("Frame rate subNyquist");
		
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