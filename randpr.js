// UNCLASSIFIED

'use strict';
/**
@class RAN
@requires stream
@requires man

Generate or learn categorical or stateless processes:

	process		#states		parms								homogeneous	?
	======================================================================
	markov 		K				TxPrs 									Y						
	bayes 		K 				eqP, net || dag						N						
	gillespie	K				states										N						
	gauss		0				[*]										N					
	wiener		0				walks										N						
	ornstein	0				theta, sigma/theta					N						

Bayes parms:

	Equilibrium probs eqP[0:K-1], net = [ [to, probs, ...], ... ] || dag = { ...}.
	
Gillespie parms:

	states = K = #states.  Transition probs TxPrs[0:K-1][0:K-1] are synthesized per 
	gillespie model
	
Markov parms:

		The K^2 (K^-K independent) transition probs TxPrs[0:K-1][0:K-1] = {0:1}, where 
		K = #states and where from-to transition probs specified to conserve prob, that is,
		must satisify sum_k TxPrs[n][k] = 1.

Gauss parms:

		values = [ ... ] ,	// pc eigen values  [unitless]
		vectors = [ [...], ...] ,	// pc eigen values	[sqrt Hz]
		ref = float,		// reference eigenvalue 
		dim = int,	// max pc dim M = observation interval T 
		mean = float // mean count in observation interval T

	These parameters generates a correlated random gaussian process and are typically 
	derived (see e.g. man) for a process with known correlation intervals M = T/Tc or 
	SNR = sqrt{ M / ( 1 + deltaC / M) }.

Wiener parms:
		walks = int // number of walks at each time step to make (stationary in first
		increments) process. 0 disables. NEEDS debug.  May need to define
		false Markov parms to init it.
		
refs:
www.statslab.cam.ac.uk/~rrw1
www.stat.yale.edu/~pollard
people.math.gatech.edu/~randall
www.stat.berkeley.edu/~pitman
www.math.dartmouth.edu/~pw
**/

var			
	// nodejs modules
	STREAM = require("stream"),		// data streams
	
	$ = require("man");   // matrix manipulators

const { EM, MVN, Copy, Each, Log } = $;
const { sqrt, floor, round, random, cos, sin, abs, PI, log, exp, min, max} = Math;

class RAN {
	
	constructor(opts, cb) {
		
		Copy({  // default configuration

			// process parameters
			wiener: 0,  // number of steps at each time step to create weiner / SSI process. 0 disables
			
			markov: null, 	// [K^2] from-to state trans probs or { states:K, index: { index: prob, ...}, ...}
				alpha: null,  // [K^2-K]/2 jump rates 
				p: null,  // [K^2-K]/2 trans probs
			
			gauss: null, // (t) => state method
			gillespie: 0, // number of states to reserve
			bayes: null, // [K] eq state probs
			
			// ensemble parameters
			
			N: 1, 		// ensemble size
			
			net: null, 	// { node: [node, ...], ... } undirected network
			dag: null, 	// { node: [node, ...], ... } directed acyclic network
			
			symbols: null, 	// state-index map = ["state1", ... "stateK"] || {state1: index, stateK: index} || K implies [0, 1, ... K]
			corrMap: null,   // map state index to correlation value [value, ... ]
			
			store: 	null,  // created by pipe()
			steps: 1, // number of process steps of size dt 
			ctmode: false, 	// true=continuous, false=discrete time mode 
			keys: null,  // event key names
			
			learn: null, 	// learner(supercb) calls back supercb(evs) || supercb(null,onend)

			emP: null, // {mu: [mean,...], sigma: [covar,...], dims: [dim, ....] } xyz-emmision probs

			filter: function (str,ev,ran) {  // filter output event ev to store/stream str
			/**
			Output event filter
				filter: function (str, ev, ran) { // event ev for stream/store str
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
			K: 0, 		// number of states (0 if stateless)
			U: null,    // [N] ensemble states (categorial) or counts (stateless)
			U0: null, 	// [N] ensemble states at t = 0
			U1: null, 	// [N] ensemble state buffer 
			UK: null,  // [N] ensemble state change accumulators (stateful K>0) or counts (stateless K=0)
			UH: null, 	// [N] ensemble holding times
			UN: null, // [N x K] ensemble counts in state ( # of times U[n] in state k )
			
			NR: null, 	// [K^2] from-to holding (mean recurrence) times
			abT: null, 	// [K'] absorption times K' <= K
			abP: null,	// [K' x K-K'] absorption probabilities K' <= K
			mleA: null, 	// [K^2] from-to state mle trans probabilities
			mleB: null, 	// {mu,sigma} observation mixing parameters
			corP: null, 	// [K^2] stat correlation probabilities
			cumP: null,	// [K^2] from-to cummulative state transition probabilities
			N0: null, 	// [K^2] from-to cummulative counts in to-state given starting from-state
			N1: null,	// [K^2] from-to state transition counts
			cumH: null, 	// [K^2] cummulative time in from-to transition
			cumN: null, 	// [K^2] cummulative number of from-to jumps
			eqP: null, 	// [K] equilibrium state probabilities 
			A: null,	// [K^2] jump rates
			
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
			emP = this.emP, // emission (aka observation) probs
			symbols = this.symbols, // state-index map
			keys = this.keys = Copy(this.keys || {}, { index:"n", state:"u", value: "v", class:"k", x:"x", y:"y", z:"z", t:"t" }); // event keys

		if ( this.alpha )  { // K-state convergent process via n = (K^2-K)/2 jump rates
			var 
				alpha = this.alpha,
				n = alpha.length,
				alpha0 = alpha.sum(),
				p = this.p = $( n, (k,p) => p[k] = alpha[k] / alpha0 );
				//Log("alpha->p", p, alpha0);
		}
		
		if ( this.p ) {   // K-state convergent process via n = (K^2-K)/2 state trans probs
			var
				p = this.p,
				n = p.length,
				K = this.K = round( (1 + sqrt(1+8*n))/2 ),
				n = 0,
				trP = this.markov = $( [K,K],  (fr,to,P) => {
					if ( to == fr ) P[fr][to] = 0;
					else
					if ( to > fr ) P[fr][to] = p[n++];
				});
			
			trP.$( (fr,to,P) => {
				if ( to == fr ) P[fr][to] = 1 - P[fr].sum();
				else
				if ( to < fr ) P[fr][to] = P[to][fr];
			});
			
			//Log("p->trP", K, trP);
		}

		this.trans = "static";	// default transition mode and refine
		
		if ( this.markov ) { // K-state process from K^2 state trans probs in K^2 - K params
			this.trans = "markov";
			var 
				trP = this.markov;

			if ( trP.constructor.name == "Object" ) {
				var
					K = this.K = trP.states,
					P = $( [K, K], $$zero);

				delete trP.states;
				for (var frKey in trP) {
					var 
						frP = trP[frKey],
						frIndex = index( frKey.split(","), [K] );

					for (var toKey in frP) {
						var toIndex = index( toKey.split(","), [K] );
						P[ parseInt(frIndex) ][ parseInt(toIndex) ] = frP[toKey];
					}
				}

				P.$$( (fr,to) =>  {
					if ( (fr==to) ) P[fr][to] = 1 - P[fr].sum();
				});
				trP = this.markov = P;
			}
			
			var
				K = this.K = trP.length,
				cumP = trP.cumP = $(K, (fr, P) => {
					P[fr] = $(K, (to, P) => {
						P[to] = trP[fr][to];
						if (to) P[to] += P[to-1];
					});
				}),
				NR = this.NR = trP.recurTimes = meanRecurTimes(trP),  // from-to mean recurrence times
				ab = trP.absorb = firstAbsorb(trP),  // first absoption times, probs, and states
				eqP = trP.eqP = $(K, (k,P) => P[k] = 1/NR[k][k]	);  // equlib state probs
		
			Log(K, trP, cumP, NR, ab, eqP);
		}
		
		else
		if ( this.bayes ) {   // bayesian network
			this.trans = "bayes";
			var 
				bayes = this.bayes,
				dag = bayes.dag, 
				net = this.net = bayes.net || dag || {},
				eqP = bayes.eqP || [0.5, 0.5],
				K = this.K = eqP.length,
				NR = this.NR = $( [K, K], ( fr, to , R ) => R[fr][to] = 1 ),
				G = bayes.G = bayes.G || $(K, (fr, G) => {  // mcmc generator G is a cummulative prob
					G[fr] = $(K, (to, G) => {
						G[to] = eqP[to];
						if (to) G[to] += G[to-1];
					});
				}), 
				V = $(N, (i, V) => { // number verticies - change this to make well-ordered or use max cardinality to make perfect numbering
					switch ("dumb") {
						case "dumb": 
					  		V[i] = i; break;
						case "well":
						case "maxcar":
					}
				}),
				A = this.A = $( [N,N], ( i, j, A ) => {	// define adjacency matrix
					if (dag)		
						deps.forEach( ( j ) => {  
							A[ j ][ i ] = 1; A[ i ][ j ] = 0;
						});
					
					else
						deps.forEach( ( j ) => {
							A[ j ][ i ] = A[ i ][ j ] = 1;
						});
				}),
				dims = net.dims = $(N, ( i, D ) => {	// dims of stores
					var deps = net[ i ] = net[ i ] || [];	// dependant nodes 
					D[ i ] = K**deps.length;
				}),
				nd = this.nd = {},	// node non-decendants
				pa = this.pa = $(N, ( i, pa ) => { //  reserve node parents
					pa[ i ] = net[ i ];
				}),
				bd = this.bd = $(N, ( i, bd ) => { // reserve boundary nodes
					bd[ i ] = new Array();
				}),
				cl = this.cl = {},	// reserve node closures
				ch = this.ch = $(N, ( i, ch ) => { // reserve children nodes
					ch[ i ] = new Array();
				}),
				alpha = net.alpha = $(N, ( i, alpha ) => {  // allocate Dirchlet priors
					var 
						alpha_i = alpha[ i ] = {};
					
					net[ i ].index( "", K, ( j ) => {	// set cond priors
						alpha_i[ j ] = $(K, ( k, alpha_ij ) => { // set all to same priors
							alpha_ij[ k ] = eqP[ k ];
						});
					});
				}),
				alpha0 = net.alpha0 = eqP.sum(),
				theta = net.theta = $(N, ( i, theta ) => { // allocate cond probs
					var 
						theta_i = theta[ i ] = {};
					
					net[ i ].index( "", K, ( j ) => {	// set cond priors
						theta_i[ j ] = $(K, ( k, theta_ij ) => { // set all to same priors
							theta_ij[ k ] = alpha[ i ][ j ][ k ] / alpha0;
						});
					});
				}),
				count = net.count = $(N, ( i, count ) => { // allocate state counters 
					var 
						count_i = count[ i ] = {};
					
					net[ i ].index( "", K, ( j ) => {	// set cond priors
						count_i[ j ] = $(K, ( k, count_ij ) => { // set all to same priors
							count_ij[ k ] = 0;
						});
					});
				});

			if ( dag )
				V.$( ( i ) => {
					net[ i ].forEach( ( j ) => {
						ch[ j ].push( i );	
					});
				});
			
			else 
				V.$( ( i ) => {
					V.$( ( j ) => {
						if ( A[ j ][ i ] ) 
							bd[ i ].push( j );
						
						else {
							var Pci = P[ i ] = {};
							Pci[ j.join(",") ] = cut( V, cl[ i ] ).join(",");
						}
					});
				});
			
			Log("bayes", N, K, dims, alpha, theta, count, G);			
		}
		
		else
		if ( this.gillespie) {	// gillespie holding time model
			this.trans = "gillespie";
			
			var
				K = this.gillespie.states,
				P = this.gillespie.P = $(K, $zero),
				NR = this.NR = $( [K, K], ( fr, to , R ) => R[fr][to] = 1 );
		}

		else
		if ( this.gauss ) {	// stateless gaussian process
			this.trans = "gauss";
		}
		
		else
		if ( this.wiener ) {	// stateless wiener process
			this.trans = "wiener";
			this.wiener.nrv = MVN( [0], [[1]] ).sample;
			Trace("wiener process experimental");
		}
		
		else
		if ( this.ornstein ) {	// stateless ornstein-ulenbeck process
			this.trans = "ornstein";
			this.ornstein.stack = $(this.steps, $zero);
		}
		
		// setup emissions
		
		if ( emP ) {	// setup emission probs
			if ( emP.dims ) {	// working on a state grid
				var 
					K = 1,
					dims = emP.dims,
					drop = dims.$( (n,Dims) => K *= Dims[n] ),
					weights = emP.weights,
					D = dims.length,
					grid = emP.grid = perms( [], dims, []),  // state grid	
					gen = emP.gen = $(K, (k,gen) => { // gauss mixing (mu,sigma) parms 
						var 
							n = 0,

							mu = $(D, (i,mu) =>
								mu[i] = grid[k][n++] + 0.5
							),

							L = $( [D,D], (i,j, L) => 	// lower trianular matrixfs with real, positive diagonal
								L[i][j] = (i <= j ) ? random() : 0
							), 

							sigma = $( [D,D], (i,j, A) => { // hermitian pos-def matrix via cholesky decomp
								var dot = 0;
								L.$( n => {
									dot += L[i][n] * L[j][n];
								});
								A[i][j] = dot * weights[i] * weights[j]
							});

						gen[k] = MVN( mu, sigma ).sample;
					});
				
				this.K = K;
			}
			
			else 
			if ( emP.mu ) {		// explicit [ (mean, covar), ....] for gaussian generators
				var
					mu = emP.mu,	// mean 
					sigma = emP.sigma || emP.cov,		// covar matricies
					K = this.K = mu.length,		// #mixes, dim(mu[k]) = D = vector dim
					parm = emP.parm = $(K, (k,p) => p[k] = {mu: mu[k], sigma: sigma[k]} ),
					gen = emP.gen = $(K, (k,g) => g[k] = MVN( parm[k].mu, parm[k].sigma ).sample );
			}
			
			else  { // derived [(mean,sigma), ....] using cholesky decomp
				var ctx = Copy( emP, {	// defaults
					decomp: "default",	// defines [ (mu,sigma), ... ] generators
					snr: 1,			// desired s/n
					cone: 0,		// randomize means on cone of angle [degs] (0=free)
					mixes: 2,		// mixes
					offcov: null,		// off-diag covars
					dim: 0,				// vector dim 
					oncov: [1,1,1]		// on-diag covars if dim=0
				});	// define gen context
				
				$({ // import these
					len: A => A._size[0],
					fiddle: (L,offcov) => {
						var k = 0,  _L = $.list(L), _offcov = $.list(offcov);
						
						return _offcov
							? $.matrix( $( L._size, (i,j, A) => A[i][j] = ( i > j ) ? _offcov[k++] : _L[i][j] ) )
							: L;
					},
					rvgen: (K,mu,sigma,cone) => {		// return K rvg generators(mu,sigma) with mean mu and covar sigma
						var 
							N = $.len(mu),	// vector dim
							rvg = {
								gen: $(K),
								parm: $(K)
							},
							xR = $.multiply( $.randRot(N), mu );
						
						const {random,cos,sin,PI} = Math;
						const {gen,parm} = rvg;
						
						Log("cone=", cone, "mixes=", K, "dim=", N);
						gen.$( n => {
							if ( cone ) {	// stay on cone of given angle from xR
								const {xG} = $(`
xR = squeeze(xR);
a = norm( xR ) * sin(theta);
N = len(xR);
V = rand(N,3);
V[:,1] = xR;
U = orthoNorm(V);
u = squeeze( U[:,2] );
v = squeeze( U[:,3] );
delta = a*( cos(phi) * u + sin(phi) * v );
xG = squeeze(xR) + delta;
`, {
	xR: $.list(xR),
	theta: 	cone*(PI/180),
	phi: 2*random() - 1
});

								parm[n] = {
									mu: $.list( xG ),
									sigma: $.list( sigma )
								};
							}
							
							else  // completely random
								parm[n] = {
									mu: $.list( n ? $.multiply( $.randRot(N), xR ) : xR ),
									sigma: $.list( sigma )
								};
							
							gen[n] = MVN( parm[n].mu, parm[n].sigma ).sample;
						});
						return rvg;
					}
				});
				
				var
					decomps = {
						default: `
N = dim ? dim : len( oncov ); 
D = diag( dim ? ones(N) : oncov ); 
L = fiddle( diag( ones(N) ), offcov );
sigma =  L * D * L'; 
mu = concat( [1], zeros(N-1) ) * snr * sqrt(sum(diag(sigma))); 
rvg = rvgen(mixes, mu, sigma, cone);
snr0 = norm(mu)/sqrt(sum(diag(sigma)));
`					},
					decomp = decomps[  ctx.decomp || "default" ] || ctx.decomp || decomps.default;
				
				const {mixes, sigma,D,L,rvg,snr0} = Copy( $( decomp, ctx), emP );
				Log(rvg);
				var K = this.K = mixes;
				emP.gen = rvg.gen;
				emP.parm = rvg.parm;
				// Log("gen snr check", snr0);
				// Log("parm", parm);
				// Log("sigma", ctx.sigma, "D", D, "L", L);
				// D = diag([4,1,9])  off=[3, -4, 5] => L = [ 1 0 0; 3 1 0; -4 5 1 ] => covar = [ 4 12 -16; 12 37 -43; -16 -43 98 ]
			}
			
			emP.obs = $(K);		// reserve observations
		}

		// define our state symbole
		
		if ( symbols )
			switch ( symbols.constructor.name ) {
				case "Object":
					K = 0;
					for (var key in symbols) K++;
					this.K = K;
					break;

				case "Array":
					for (var syms = {}, k=0, K=this.K=symbols.length; k<K; k++) syms[ symbols[k] ] = k;
					symbols = this.symbols = syms;
					break;

				default:
					for (var syms = {}, k=0, K=this.K=symbols; k<K; k++) syms[ k ] = k;
					symbols = this.symbols = syms;
			}	
		
		else {
			for (var syms = {}, k=0, K=this.K; k<K; k++) syms[ k ] = k;
			symbols = this.symbols = syms;
		}
		
		// create zero-mean correlation map
		
		var map = this.corrMap = new Array(K);
		
		if (K) 
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

		Log({
			keys: keys,
			states: K, 
			syms: symbols, 
			xmap: map,
			emP: emP
		});	

		// allocate ensemble vars and state counters
		
		var 
			U1 = this.U1 = $(N),
			UK = this.UK = $(N, $zero),
			N1 = this.N1 = $( [K,K], $$zero),	
			mleA = this.mleA = $( [K,K], $$zero), 
			cumH = this.cumH = $( [K,K], $$zero),
			cumN = this.cumN = $( [K,K], $$zero),
			mleR = this.mleR = $( [K,K] ),
			err = this.err = 1,
			corP = this.corP = $( [K,K] ),
			emP = this.emP,
			eqP = this.eqP = $(K, $zero),
			p = 1/K,
			Np = p * N,
			N0 = this.N0 = $( [K,K],  (fr,to,N0) => N0[fr][to] = (fr == to) ? Np : 0 ),
			UH = this.UH = $(N),
			U = this.U = $(N),
			U0 = this.U0 = $(N),
			ctmode = this.ctmode,
			UN = this.UN = $( [N, K], $$zero);
		
		this.t = this.s = this.samples = 0;  // initialize process counters

		// initialize ensemble
		
		if ( this.learn )   // in learning mode
			U.$( n => UH[n] = U0[n] = U[n] = 0 );
		
		else
		if (K) { // initialize K-state process
			if (this.markov) {	// init markov process
				if ( K == 2) { // special 2-state process
					var R01=NR[0][1], R10=NR[1][0];

					U.$( n => {
						if ( n < Np ) {
							var fr = U0[n] = U[n] = 1;
							UH[n] = NR[fr][fr] = ctmode ? expdev(-1/A[fr][fr]) : 0;
						}

						else {
							var fr = U0[n] = U[n] = 0;
							UH[n] = NR[fr][fr] = ctmode ? expdev(-1/A[fr][fr]) : 0;
						}
					});
				}

				else // general K-state process
					U.$( n => {
						var fr = floor(random() * K);
						U0[ n ] = U[n] = fr; 
						UH[ n ] = NR[fr][fr] = ctmode ? expdev(-1/A[fr][fr]) : 0;
						UN[ n ][ fr ] = 1;
					}); 
			}
			
			else   // initialize stateful process
				U.$( n => { UH[n] = 0; U0[n] = U[n] = 0; } );
		}
		
		else	// init stateless process
			U.$( n => { UH[n] = 0; U0[n] = U[n] = 0; } );
		
		//Log("UH", UH);

		this.gamma = $(this.steps, $zero);
		
		if (cb) cb(null);
	}
	
	statCorr( ) {  // statistical correlation function
		var 
			K = this.K, map = this.corrMap, cor = 0, corP = this.corP, p, N0 = this.N0, N = this.N, samples = this.samples;

		if (samples)
			map.$( fr => {
				map.$( to => {
					p = corP[fr][to] = N0[fr][to] / samples;
					cor += map[fr] * map[to] * p;
				});
			});
		
		else
			cor = 1;
		
		this.samples += N;

		return cor ; 
	}

	eqProbs( ) { // equlib probs
		var
			N = this.N,
			T = this.t,
			eqP = this.eqP,
			UN = this.UN;
		
		eqP.$( (k,P) => {
			var sum =0;
			UN.$( (n,N) => sum += N[k] );
			P[k] = sum / N / T;
		});
	}
	
	condProbs( ) {	// estimate network conditional probs from Directlet priors
		var 
			K = this.K,
			U = this.U,
			UN = this.UN,
			net = this.net,		
			dims = net.dims,
			alpha = net.alpha,
			theta = net.theta,
			count = net.count;

		U.$( i => {
			var 
				j = net[ i ].index( U ),
				counts = count[ i ][ j ],
				thetas = theta[ i ] [ j ],
				alphas = alpha[ i ] [ j ],
				Ucounts = UN[ i ];

			Ucounts.$( k => {
				counts[ k ] += Ucounts[ k ];
			});

			var
				count0 = counts.sum(),
				alpha0 = alphas.sum();

			Ucounts.$( k => {
				thetas[ k ] = ( counts[ k ] + alphas[ k ] ) / ( count0 + alpha0 );
			});
		});
	}

	holdTimes( ) {	// mle hold times
		var
			cumH = this.cumH,
			cumN = this.cumN,
			mleR = this.mleR;
		
		mleR.$( fr => {   // estimate jump rates using cummulative UH[fr][to] and N[fr][to] jump times and counts
			mleR[fr].$( to => {
				mleR[fr][to] = (fr == to) ? 0 : cumH[fr][to] / cumN[fr][to];
			});
		});
	}
	
	transProbs( ) {	// mle transition probs
		var
			N1 = this.N1,
			mleA = this.mleA;
		
		N1.$( fr => {  // transition prob mleA using the 1-step state transition counts
			N1[fr].sum( sum => {
				mleA[fr].$( to => {
					mleA[fr][to] = N1[fr][to] / sum;
					//Log(fr,to,Nfr[to], sum);
				});
			});
		});
	}		
	
	step (evs, cb) {  // advance process forward one step (with events evs if in learning mode)
		
		function draw( P ) { // draw random state with cumulative prob P
			var to = 0, K = P.length;

			for (var u = random(); to < K && P[to] <= u; to++) ;
			return (to == K) ? to-1 : to;
		}

		var 
			ran = this,
			UH = this.UH, NR = this.NR, mleA = this.mleA,
			
			UK = this.UK, U1 = this.U1, U = this.U, U0 = this.U0, UN = this.UN,
			N1 = this.N1, N0 = this.N0,
			cumH = this.cumH, cumN = this.cumN, A = this.A, 
			
			symbols = this.symbols, keys = this.keys, emP = this.emP, 
			K = this.K, t = this.t, N = this.N, s=this.s, dt = this.dt,
			
			transitions = {
				// homogeneous stateful process

				static: function ( t, u ) { 
					return u;
				},
				
				markov: function ( t, u ) {  // toState via trans probs
					var 
						markov = ran.markov,
						cumP = markov.cumP,
						fr = u,
						to = draw( cumP[fr] );

					return to;
				},

				// inhomogeneous stateful process

				gillespie: function ( t, u ) {  // toState via trans probs computed on holding times
					var 
						gillespie = ran.gillespie,
						cumP = gillespie.P,
						fr = u,
						R0 = NR[fr], 
						K = cumP.length;

					cumP.$( (to, P) => {
						P[to] = (fr==to) ? 0 : NR[to] / R0;
						if (to) P[to] += P[to-1];
					});

					var P0 = cumP[K-1];
					cumP.$( (to,P) => P[to] /= P0 );

					return draw( cumP[fr] );
				},

				bayes: function ( t, u ) {  // toState using metropolis-hastings (or mcmc) with generator G
					var 
						bayes = ran.bayes,
						P = bayes.eqP,
						G = bayes.G,
						fr = u,
						toG = G[fr], 
						to = draw(toG),
						frG = G[to],
						Ap = P[to] / P[fr],
						Ag = frG[to] / toG[fr],
						accept = min(1 , Ap * Ag), // acceptance criteria
						u = random();

					//if (t == 0) Log(">>>>>>>>>>>Ginit", t, G, P, accept);
					return (u <= accept) ?  to : fr;
				},

				// stateless process

				gauss: function ( t, u ) {
					var
						gauss = ran.gauss,
						vals = gauss.values,	// pc eigen values  [unitless]
						vecs = gauss.vectors,	// pc eigen values	[sqrt Hz]
						ref = gauss.ref,	// ref eigenvalue
						N = vals.length, // number of eigenvalues being used
						dim = gauss.dim,	// max pc dim = observation interval
						mean = gauss.mean;  // mean events over sample time

					if ( t >= dim ) // KL expansion exhausted
						return mean;	// could return negbin dev but good to know it is no longer updating

					else {	// use KL expansion to generate count sample from arrivial rate process 
						var
							B = $.matrix( $(N, (n,B) => {  // generate KL coefficients [events]
								var 
									Bmod = sqrt( expdev( mean * vals[n] / ref ) ),  
									Barg = random() * PI;

								//if (t == 0)  Log( t , n , N, vals[n] / ref , mean, Bmod);
								B[n] = $.complex( Bmod * cos(Barg), Bmod * sin(Barg) );  
							}) ),
							V = $.matrix( $(N, (n,V) => {  // get KL eigen vector at time t [sqrt Hz]
								V[n] = vecs[n][t];  // t index = 0:N samples vectors at t = -T/2 : T/2
							}) ),
							A = $.dot( B, V),  // complex analytic function representing event rate process [sqrt Hz]
							lambda = $.abs(A)**2,  // event rate process [Hz]
							k = lambda * dt;  // [events]

						//if (t == 0 ) Log("gauss", t, N, dim, mean, k, ref, A);
						return k;
					}
				},

				wiener: function ( t, u ) {
					var 
						wiener = ran.wiener,
						M = wiener.walks,
						nrv = wiener.nrv;

					for (var j=1, walks=floor(M*t); j<=walks; j++) u += nrv()[0];

					return u / sqrt(M);
				},

				ornstein: function ( t, u ) {
					var 
						ornstein = ran.ornstein,
						theta = ornstein.theta,
						a = ornstein.a,  // a = sigma / sqrt(2*theta)
						stack = ornstein.stack,
						Et = exp(-theta*t),
						Et2 = exp(2*theta*t),
						Wt = stack[ floor(Et2 - 1) ] || 0;

					stack.push( Wt );

					return a * Et * Wt;	
				}				
			},

			trans = transitions[ ran.trans ];
				
		U.$( n => U1[n] = U[n] );  // hold current states/values
		
		if (evs) { // in learning mode 
			/*
			if ( !t ) { // initialize states at t=0
				ran.gamma[0] = 1;  // force
				evs.forEach( ev => {
					var n = ev[keys.index];		// ensemble index
					U[ n ] = symbols[ev[keys.state]] || 0;  // state (0 if hidden)
				});
				U.$( n => {
					U0[n] = U1[n] = U[n];		// initialize states
					UK[n] = -1;  // force initial jump counter to 0
				}); 
			} */

			// latch (assume time-ordered) events to ensemble
			
			t = this.t = evs[0].t;	
			
			if (K) // categorical process so latch states
				evs.forEach(ev => {	// set states (if supervised) or symbols[0] (if hidden)
					U[ ev[keys.index] || 0 ] = symbols[ ev[keys.value] || 0 ];
				});
			
			else // stateless process so latch values
				evs.forEach(ev => {  // set values (if supervised) or 0 (if hidden)
					U[ ev[keys.index] || 0 ] = ev[keys.value] || 0;
				});				
		}

		else  { // in generative mode
			U.$( n => {
				U[ n ] = trans( t , U[n] );
			});
		}
		
		// Log(t, U.join(" "));
		
		// update process counters
		
		if (K)  { // categorical process
			this.gamma[s] = this.statCorr();		

			U.$( n => {  // adjust jump counters
				var
					frState = U1[n],
					toState = U[n];

				if ( frState != toState) { // jump if state changed
					var
						held = t - UH[n],	// initially 0 and remains 0 in discrete-time mode
						hold = this.ctmode ? expdev( 1/A[frState][toState] ) : 0 ;  // draw expected holding time

					cumH[frState][toState] += held; // cummulative holding time in from-to jump
					cumN[frState][toState] ++;  // cummulative number of from-to jumps
					
					NR[frState][frState] = hold;  // update expected holding time 
					UH[ n ] = t + hold;    // advance to next jump time (hold is 0 in discrete time mode)
					UK[ n ]++; 		// increment jump counter

					ran.onJump(n, toState, hold ); 	// callback with jump info
				}
			});	
			
			U.$( n => {   // adjust state counters
				var k = U[n]; 				// state
				N0[ U0[n] ][ k ]++; 	// initial-to counters for computing ensemble correlations
				N1[ U1[n] ][ k ]++;		// from-to counters for computing trans probs
				UN[ n ] [ k ]++; 		// # times U[n] in state k; for computing cond probs
			});
			
			if ( emP ) 
				U.$( n => emP.obs[n] = emP.gen[ U[n] ]() );
		}
	
		else  // stateless process
			U.$( n => {   // adjust counters
				UK[ n ] += U[ n ];
			});
			
		//Log( (t<10) ? "0"+t : t, U.join(""));
		//if (t<50) Log( t<10 ? "0"+t : t,U,UK);		
		//if (t<5) Log(t,N0);
		
		this.onStep();
		this.t += this.dt;
		this.s++;
	}
	
	start ( ) {	  // start process in learning (reverse) or generative (forward) mode
		var 
			ran = this,
			U = this.U,
			interpolate = false,
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
							store: ran.store,  // output event store
							T: ran.steps,		// observation time
							F: ran.F,	// event count frequencies
							f: mass( ran.F, N, false),	// event count probabilities
							Nevs: counts( ran.F ), 	// total number of events
							J: ran.UK,		// ensemble counts
							N: ran.N		// ensemble size
						});
				}

				if ( batch )
					if ( ran.s % batch == 1 ) ran.onBatch();
			});
		
		else { // generative mode
			//Log("start gen", ran.steps, ran.N);
			
			while (ran.s < ran.steps) {  // advance process to end
				ran.step(null);
				
				if ( batch )
					if ( ran.s % batch == 1 ) ran.onBatch();
			}
			
			ran.onEnd();
		}
		
	}
	
	corrTime ( ) {  // return correlation time computed as area under normalized auto correlation function
		
		if ( this.K ) { // categorical process
			var Tc = 0;
			for (var t=0, T = this.t; t<T; t++) Tc += abs(this.gamma[t]) * (1 - t/T);

			Tc *= this.dt / this.gamma[0] / 2;
		}
		
		else
			var Tc = 0;
		
		Log(">>>>>>>>>Tc=", Tc);
		return Tc;
	}
	
	countFreqs ( ) {  // return count frequencies across the ensemble

		var
			UK = this.UK,  // ensemble counters
			K = this.K,	// number of states (0 = stateless process)
			Kmax = K ? UK.max() : floor( UK.max() ),
			F = this.F = $( 1 + Kmax , $zero );  

		UK.$( n => F[ K ? UK[n] : floor( UK[n] ) ]++ ); 
		return F;
	}
	
	record (at, ev) {  // record event ev labeled at to store or stream
		ev.t = this.t;
		ev.at = at;
		this.filter(this.store, ev, this);
	}
	
	// State methods
	
	onBatch () {    // record MLE jump rates and trans probs
		var 
			ran = this,
			K = this.K,
			t = this.t,
			s = this.s,
			N = this.N,
			net = this.net,
			F = this.countFreqs();

		if (K) {
			if ( this.net ) this.condProbs( );
			this.eqProbs( );
			this.transProbs( );
			this.holdTimes( );
		}
		
		else
		if ( this.markov ) {	// relative error between mle and actual trans probs
			var 
				trP = this.markov,
				ref = 0,
				err = this.err = trP   
					? abs( mleA[ref][ref] - trP[ref][ref] ) / trP[ref][ref]
					: 0;
		}
		
		else
			var err = 0;
		
		Log("batch", t, F.length, this.UK.avg().toFixed(4), net ? net.theta : null );
		//Log("batch", t, F.length, this.UK.avg().toFixed(4), F.join(" "));
		
		this.record("batch", {
			count_freq: F,
			count_prob: $( F.length, n => F[n]/N ),
			rel_error: err,
			mle_tr_probs: this.mleA,
			mle_hold_time: this.mleR,
			eq_probs: K ? this.eqP : null,
			eq_cond_probs: net ? net.theta : null,
			stat_corr: K ? this.gamma[ s-1 ] : 0
		});	
	}

	onError ( msg ) {	// record process error condition
		Log(msg);
		this.record("error", { 
			error: msg
		});
	}
	
	onJump (index,state,hold) {  // record process jump info		
		this.record("jump", {
			index: index, state:state, hold:hold
		});
	}
	
	onStep () {		// record process step info
		this.record("step", {
			gamma:this.gamma[this.s],
			walk: this.wiener ? this.U : []
		});
	}

	onConfig() {  // record process config info
		this.record("config", {
			states: this.K,
			ensemble_size: this.N,		
			sample_time: this.dt,
			nyquist: 1/this.dt,
			cum_tr_probs: this.cumP,
			
			markov_tr_probs: this.markov,
			//trans_mode: this.transMode,
			
			mean_recur_times: this.NR,
			eq_probs: this.eqP,
			mixes: this.emP ? this.emP.parm : null,
			snr: this.emP ? this.emP.snr0 : 0,
			run_steps: this.steps,
			absorb_times: this.ab
		});
	}
	
	onEnd() {  // record process termination info
		
		//Log("onend", this.obs.length);
		
		var 
			ran = this,
			batch = this.batch,
			T = this.steps,
			Tc = this.corrTime(),
			Kbar = this.UK.avg(),
			M = T / Tc,
			delta = Kbar / M,
			F = this.countFreqs(),
			K = this.K,
			mleB = this.mleB = this.emP ? EM( this.emP.obs, K) : null;

		//Log("onend UK", UK);
		
		this.record("end", {  // record supervised stats
			stats: {
				mle_holding_times: ran.mleR,
				rel_error: ran.err,
				count_freq: F,
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

		if (false)
			if (evs = this.store)
				evs.forEach( ev => {
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

function counts(F) {
	var Nevs = 0;
	F.$( k => Nevs += k * F[k] );
	return Nevs;	
}

function mass(F, N, interpolate) {
	return $( F.length, (k, f) => {    // count probabilities
		if (interpolate)  {  // usually pointless
			if ( F[k] ) 
				f[k] = F[k] / N;

			else
			if ( k ) {
				N += F[k-1];
				f[k] = F[k-1] / N;
			}

			else
				f[k] = 0;
		}
		
		else
			f[k] = F[k] / N;
	});
}

function cumulative( P ) {  // replace P with its cumulative
	switch (0) {
		case 0:
			P.$( (k) => {
				if (k) P[k] += P[k-1];
			});
			break;
			
		case 1:
			var cum = 0;
			P.$( (k) => {
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
	A.$( (k) => A[k][k] = - A[k].sum() );
	return A;
}

function balanceProbs(P) {  // enforce global balance on probs
	P.$( (k) => {
		P[k][k] = 1 - P[k].sum() 
	});
	return P;
}			

function Trace(msg) {
	"rand".trace(msg);
}

function firstAbsorb(P) {  //< compute first absorption times, probs and states
	var 
		K = P.length,
		kAb = [],
		kTr = [],
		x = P.$( (k) => {
			if ( P[k][k] == 1 ) 
				kAb.push(k+1);
			else
				kTr.push(k+1);
		}),
		ctx = {
			P: $.matrix(P),
			K: K,
			kAb: $.matrix(kAb),
			kTr: $.matrix(kTr),
			nAb: kAb.length,
			nTr: kTr.length,
			abT: $.matrix([]),
			abP: $.matrix([])
		};
	
	//Log("ab ctx", JSON.stringify(ctx));
	if ( ctx.nAb && ctx.nTr )
		$("Q = P[kTr,kTr]; NR = P[kTr,kAb]; N = inv( eye(nTr,nTr) - Q ); abT = N*ones(nTr,1); abP = N*NR;", ctx);
		
	return {
		times: ctx.abT,
		probs: ctx.abP,
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
			P: $.matrix(P),
			K: P.length
		},
		K = ctx.K;

	if ( K > 1) {
		$("k=2:K; P0=P[1,1]; Pl=P[k,1]; Pu=P[1,k]; Pk=P[k,k]; A = Pk - eye(K-1); Adet = abs(det(A)); ", ctx);

		Log({"MRT det": ctx.Adet});

		if ( ctx.Adet < 1e-3 ) {
			Log("Proposed process is not ergodic, thus no unique eq prob exist.", ctx.P);
			return $( [K,K],  $$zero );
		}

		else {
			$("wk= -Pu*inv(A);", ctx);

			//Log("man meanRecurTimes", ctx);
			ctx.w = $.matrix([1].concat(ctx.wk[0]));

			$("w = w / sum(w); w = [w]; Z = inv( eye(K) - P + w[ ones(K) , 1:K] ); H = zeros(K,K); ", ctx);

			var 
				H = ctx.H,
				Z = ctx.Z,
				w = ctx.w[0];

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
	
	grid.$( n => {
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
	function index(key, K, cb) {
		if (cb)
			if (key.length == this.length) 
				cb( key );

			else
				for (var k=0; k<K; k++)
					this.index( key+k, K, cb );
		
		else {
			var 
				vars = key,
				key = "";
			
			for (var n=0, N=this.length; n<N; n++) key += vars[ this[n] ];

			return key;
		}
	},
	
	function sample(delta) {
		var 
			A = this,
			k = 0,
			rtn = $( floor(A.length/delta), (n,NR) => NR[n] = A[k += delta] );
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
		A.$( (k) => {
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
].Extend(Array);

switch ( process.argv[2] ) {   //< unit tests
	case "?":
		Log("unit test with 'node randpr.js [R1 || R2 || ...]'");
		break;
		
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
			[[1,0,0,0,0],[0.5,0,0.5,0,0],[0,0.5,0,0.5,0],[0,0,0.5,0,0.5],[0,0,0,0,1]] // 2 absorbing states
/*
{ times: [ [ 3 ], [ 3.9999999999999996 ], [ 2.9999999999999996 ] ],
  probs: 
   [ [ 0.75, 0.24999999999999997 ],
     [ 0.49999999999999994, 0.49999999999999994 ],
     [ 0.24999999999999997, 0.7499999999999999 ] ],
  states: [ 1, 5 ] }
*/
			//[[0.5,0.25,0.25],[0.5,0,0.5],[0.25,0.25,0.5]]  // no absorbing stats
/*
{ times: [], probs: [], states: [] }
*/
		));
		break;
		
	case "R2.1":  // config methods
		var ran = new RAN({
			p: [.4],
/*
p->trP 2 [ [ 0.6, 0.4 ], [ 0.4, 0.6 ], rows: 2, columns: 2 ]
R> { 'MRT det': 0.4 }
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
  xmap: [ 1, -1 ] }
*/

			//markov: [[0.1, 0.9], [0.1, 0.9]]
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
  xmap: [ 1, -1 ],
  trProbs: [ [ 0.1, 0.9 ], [ 0.1, 0.9 ] ] }
  */
			//markov: { states: 3, 0: {1: 0.8, 2: 0.1}, 1: {0: 0.1} }
/*
R> { 'MRT det': 0 }
R> Proposed process is not ergodic, thus no unique eq prob exist. Matrix {
  _data: 
   [ [ 0.09999999999999998, 0.8, 0.1 ],
     [ 0.1, 0.9, 0 ],
     [ 0, 0, 1 ],
     rows: 3,
     columns: 3 ],
  _size: [ 3, 3 ],
  _datatype: undefined }
R> { keys: 
   { index: 'n',
     state: 'u',
     class: 'k',
     x: 'x',
     y: 'y',
     z: 'z',
     t: 't' },
  states: 3,
  syms: { '0': 0, '1': 1, '2': 2 },
  xmap: [ 0, 1, -1 ] }
  */
		});
		break;

	case "R2.3":  // config methods
		var ran = new RAN({
			emP: {
				dims: [3,3],
				weights: [1,1]
			},
			markov: { states: 9, 0: {1: 0.8, 2: 0.1}, 1: {0: 0.1}, "0,1": { "1,0": .4} }
		});
/*
markov [ [ 0.09999999999999998, 0.8, 0.1, 0, 0, 0, 0, 0, 0 ],
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
			markov: "random"
		});
		break;
		
	case "R3":  // sync pipe with various textbook examples, custom filtering with supervised learning validation
		var ran = new RAN({
			// these have same eqprs [.5, .5] (symmetry -> detailed balance --> eqP[k] = 1/K  eqpr)
			//markov: [[.6, .4],[.4, .6]],
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
MLEs { holdTimes: '[[0,2.4961510304589924],[2.4806757386729417,0]]',
  emProbs: 'null',
  trProbs: '[[0.601384811290722,0.39861518870927803],[0.40071697978935383,0.5992830202106462]]' }
*/
			//markov: [[0.83177, 0.16822], [0.17152, 0.82848]],
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
MLEs { holdTimes: '[[0,5.886354079058032],[5.810975609756097,0]]',
  emProbs: 'null',
  trProbs: '[[0.832333007235839,0.167666992764161],[0.169801047868115,0.830198952131885]]' }
  */
			//markov: [[.5, .5], [.5, .5]],
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
MLEs { holdTimes: '[[0,1.9922586601800394],[1.9936797628135707,0]]',
  emProbs: 'null',
  trProbs: '[[0.5003326916081868,0.4996673083918131],[0.49942810529955756,0.5005718947004424]]' }
  */
			//markov: [[0.1, 0.9], [0.9, 0.1]],

			// textbook exs
			markov: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3
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
MLEs { holdTimes: '[[0,1.1008115157955753],[9.660927669121593,0]]',
  emProbs: 'null',
  trProbs: '[[0.10186550976138829,0.8981344902386117],[0.10127438873795999,0.89872561126204]]' }
*/
			
			//markov: [[1/2, 1/3, 1/6], [3/4, 0, 1/4], [0,1,0]],  // pg142 ex2  eqpr [.5, .333, .1666]
			//markov: [[1,0,0], [1/4, 1/2, 1/4], [0,0,1]],  // pg143 ex8  no eqprs

			// these have different eqprs
			//markov: [[0.9,0.1],[0.1,0.9]],
			//markov: [[0.1, 0.9], [0.1, 0.9]],  // bernoulli scheme has identical rows
			//markov: [[0.1, 0.9], [0.3, 0.7]],
			//markov: [[0.1, 0.9], [0.4, 0.6]],

			// textbook exs 
			//markov: [[0,1],[1,0]],  // pg433 ex16  regular (all states reachable) absorbing/non on even/odd steps non-regular non-absorbing but ergodic so --> eqpr [.5, .5]
			//markov: [[0.5,0.25,0.25],[0.5,0,0.5],[0.25,0.25,0.5]],  // pg406 ex1  regular (after 2 steps) thus ergodic so eqpr [.4, .2, .4]
			//markov: [[0,1,0,0,0], [0.25,0,0.75,0,0], [0,0.5,0,0.5,0], [0,0,0.75,0,0.25], [0,0,0,1,0]],  // pg433 ex17  non-absorbing non-regular but ergodic so eqpr [.0625, .25, .375]
			//markov: [[1,0,0,0,0],[0.5,0,0.5,0,0],[0,0.5,0,0.5,0],[0,0,0.5,0,0.5],[0,0,0,0,1]],    // 2 absorbing states; non-ergodic so 3 eqpr = [.75 ... .25], [.5 ... .5], [.25 ...  .75]

			//markov: [[1-.2, .1, .1], [0, 1-.1, .1], [.1, .1, 1-.2]],
			//markov: [[1-.2, .1, .1], [0.4, 1-.5, .1], [.1, .1, 1-.2]],
			//markov: [[1-.6, .2, .2,.2], [.1, 1-.3, .1,.1], [.1, .1, 1-.4,.2],[.1,.1,1-.8,.6]],  // non-ergodic
			
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
MLEs { holdTimes: '[[0,1.0968340824701974],[9.798737174427782,0]]',
  emProbs: '[
  	{"weight":0.5444315148931079,"mu":[1.0139285640388083],"sigma":[[1.020168448272447]],"_gaussian":
  		{"sigma":[[1.020168448272447]],"mu":[1.0139285640388083],"k":1,"_sinv":[[0.9802302763758278]],"_coeff":0.3949791055911031}},
		
	{"weight":0.4555684851068921,"mu":[1.071684787418599],"sigma":[[2.070246038073367]],"_gaussian":
		{"sigma":[[2.070246038073367]],"mu":[1.071684787418599],"k":1,"_sinv":[[0.48303437447011366]],"_coeff":0.2772675754216858}}]',
  trProbs: '[[0.09848484848484848,0.9015151515151515],[0.09996252391565909,0.9000374760843409]]' }
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
							holdTimes: JSON.stringify(H),
							emProbs: JSON.stringify(B),
							trProbs: JSON.stringify(A)
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
		
	case "R3.1":  // gen process for R3.2 with async pipe to stdout
		var ran = new RAN({

			markov: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3

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
		
	case "R3.2":  // gen process for R3.3 using async pipe to stdout
		var ran = new RAN({

			markov: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3

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
		
	case "R3.3":  // supervised learning with R3.2 evs using sync pipe to store
		var 
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
					evs.$( "group", (evs) => {
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

				markov: {},  
				//keys: {state:"u", index: "n"},
				K: 2,  // assume 2-state process
				N: 50  // assume 50 members in ensemble
			});
	
		ran.pipe( function (store) {
			Log(store);
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
