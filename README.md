/**
@class RANDPR
	[SourceForge](https://sourceforge.net) 
	[github](https://github.com/acmesds/debe.git) 
	[geointapps](https://git.geointapps.org/acmesds/debe)
	[gitlab](https://gitlab.weat.nga.ic.gov/acmesds/debe.git)
	
# RANDPR

In its markov generation mode, RANDPR will generate a first-order, K-state Markov process given 
its first-order, KxK transition probabilites.  RANDPR also computes equilibrium probabilities (if they 
exist), mean recurrence times, time to first absorption, absorption probabilites and the statistical 
auto-covariance function.  RANDPR can create ergodic, regular, absorptive processes, Weiner 
(stationary in first increments), and wide-sense stationary processes.

In its markov learning mode, RANDPR produces supervised and unsupervised estimates:
MLEs of the underlying transition probabilities, number of coherence intervals (and related SNR),
and the underlying intensity profile.  (MLE for the Weiner process, e.g. first time to exit, have not 
yet been implemented).

Both discrete- and continious-time models are supported in either forward or reverse mode.  

RANDPR can be customized with onStep(), onBatch(), onEnd(), and filter() methods to 
sample metrics in both forward and reverse modes.  Both modes can make use of piped 
streams to minimize memory usage.

RANDPR also supports ornstein, bayes, gillespie, and weiner processes.

## Installation

Clone [RANDPR random process](https://github.com/acmesds/randpr) into your PROJECT/randpr folder.  
Clone [ENUM basic enumerators](https://github.com/acmesds/enum) into your PROJECT/enum folder.  
Clone [MAN matrix manipulator](https://github.com/acmesds/jslab) into your PROJECT/man folder.  

## Usage

Require the module:

	var RAN = require("randpr");
	
then create a new instance:

	var ran = new RAN({
		key: value, 						// set key
		"key.key": value, 					// indexed set
		"key.key.": value					// indexed append
	}, function (err) {
		console.log( err ? "something evil is lurking" : "look mom - Im running!");
	});

where [its configuration keys](/shares/prm/randpr/index.html) follow 
the [ENUM deep copy conventions](https://github.com/acmesds/enum).

The instance can then be piped:

	ran.pipe()
	
to genenerate a process (if configured in the forward/generate mode), or learn process parameters (if 
configured in the reverse/learning mode).

The following examples are from RANDPR's unit tester:

	node randpr.js [R1 || R2 || ... ]

### R2.1 - config methods
	var ran = new RAN({
		p: [.4],
		//markov: [[0.1, 0.9], [0.1, 0.9]]
		//markov: { states: 3, 0: {1: 0.8, 2: 0.1}, 1: {0: 0.1} }
	});

### R2.3 - config methods
	var ran = new RAN({
		emP: {
			dims: [3,3],
			weights: [1,1]
		},
		markov: { states: 9, 0: {1: 0.8, 2: 0.1}, 1: {0: 0.1}, "0,1": { "1,0": .4} }
	});
		
### R2.4 - config methods
	var ran = new RAN({
		emP: {
			dims: [2,2,2],
			weights: [1,1,1]
		},
		markov: "random"
	});
		
### R3 sync pipe with various textbook examples, custom filtering with supervised learning validation
	
	var ran = new RAN({
		// these have same eqprs [.5, .5] (symmetry -> detailed balance --> eqP[k] = 1/K  eqpr)
		//markov: [[.6, .4],[.4, .6]],

		//markov: [[0.83177, 0.16822], [0.17152, 0.82848]],

		//markov: [[.5, .5], [.5, .5]],
		//markov: [[0.1, 0.9], [0.9, 0.1]],

		// textbook exs
		markov: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3

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
			Log(store);
		});
		
### R3.1 - gen process for R3.2 with async pipe to stdout
	
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
		
### R3.2 - gen process for R3.3 using async pipe to stdout

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
	
	ran.pipe(process.stdout);  // stdout evs used in R3.3
		
### R3.3 - supervised learning with R3.2 evs using sync pipe to store
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
				evs.$( true, (evs) => {
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
	
## Contributing

See our [issues](https://totem.west.ile.nga.ic.gov/issues.view), [milestones](https://totem.west.ile.nga.ic.gov/milestones.view), [s/w requirements](https://totem.west.ile.nga.ic.gov/swreqts.view),
and [h/w requirements](https://totem.west.ile.nga.ic.gov/hwreqts.view).
	
## License

[MIT](LICENSE)

*/