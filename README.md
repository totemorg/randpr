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
	
to genenerate a process (if configured in the forward/generate mode), or learn the process parameters (if 
configured in the reverse/learning).

### Simple config
	var ran = new RAN({
		trP: [[0.1, 0.9], [0.1, 0.9]],  
	});
		
### Sparse config
	var ran = new RAN({
		K: 3,
		trP: { 0: {1: 0.9}, 1: {0: 0.1} }
	});
	
### Sparse config with 2D emission parameters
	var ran = new RAN({
		K: 9,
		obs: {
			dims: [3,3],
			weights: [1,1]
		},
		trP: { 0: {1: 0.9}, 1: {0: 0.1}, "0,1": { "1,0": .4} }
	});

### Random transition probs with 3D emission parameters
	var ran = new RAN({
		K: 8,
		obs: {
			dims: [2,2,2],
			weights: [1,1,1]
		},
		trP: "random"
	});
	
### forward mode with textbook examples, custom filtering, supervised and unsupervised learning
	var ran = new RAN({
		// these have same eqprs [.5, .5] (symmetry -> detailed balance --> pi[k] = 1/K  eqpr)
		//trP: [[.6, .4],[.4, .6]],
		//trP: [[0.83177, 0.16822], [0.17152, 0.82848]],
		//trP: [[.5, .5], [.5, .5]],
		//trP: [[0.1, 0.9], [0.9, 0.1]],

		// textbook exs
		trP: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3
		//trP: [[1/2, 1/3, 1/6], [3/4, 0, 1/4], [0,1,0]],  // pg142 ex2  eqpr [.5, .333, .1666]
		//trP: [[1,0,0], [1/4, 1/2, 1/4], [0,0,1]],  // pg143 ex8  no eqprs

		// these have different eqprs
		//trP: [[0.9,0.1],[0.1,0.9]],
		//trP: [[0.1, 0.9], [0.1, 0.9]],  // bernoulli scheme has identical rows
		//trP: [[0.1, 0.9], [0.3, 0.7]],
		//trP: [[0.1, 0.9], [0.4, 0.6]],

		// textbook exs 
		//trP: [[0,1],[1,0]],  // pg433 ex16  regular (all states reachable) absorbing/non on even/odd steps non-regular non-absorbing but ergodic so --> eqpr [.5, .5]
		//trP: [[0.5,0.25,0.25],[0.5,0,0.5],[0.25,0.25,0.5]],  // pg406 ex1  regular (after 2 steps) thus ergodic so eqpr [.4, .2, .4]
		//trP: [[0,1,0,0,0], [0.25,0,0.75,0,0], [0,0.5,0,0.5,0], [0,0,0.75,0,0.25], [0,0,0,1,0]],  // pg433 ex17  non-absorbing non-regular but ergodic so eqpr [.0625, .25, .375]
		//trP: [[1,0,0,0,0],[0.5,0,0.5,0,0],[0,0.5,0,0.5,0],[0,0,0.5,0,0.5],[0,0,0,0,1]],    // 2 absorbing states; non-ergodic so 3 eqpr = [.75 ... .25], [.5 ... .5], [.25 ...  .75]

		//symbols: [-1,1],

		//trP: [[1-.2, .1, .1], [0, 1-.1, .1], [.1, .1, 1-.2]],
		//trP: [[1-.2, .1, .1], [0.4, 1-.5, .1], [.1, .1, 1-.2]],
		//symbols: [-1,0,1],
		//trP: [[1-.6, .2, .2,.2], [.1, 1-.3, .1,.1], [.1, .1, 1-.4,.2],[.1,.1,1-.8,.6]],  // non-ergodic

		solve:  {
			batch: 50,  // supervised learning
			compress: true,
			interpolate: false,
			lma: [50]  // unsuervised learning
			// bfs: [5,200,5]
			// lfa: [50]
		},

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
					//Log(ev);
					str.push(ev);
					break;
			}
		},

		N: 1000,
		steps: 200
	}, function (ran) {
		ran.pipe( function (store) {
			Log(store);
		});
	});
		
### Async pipes
	var ran = new RAN({

		trP: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3

		solve:  {
			batch: 50,
			compress: true,
			interpolate: false,
			lma: [50]
			// bfs: [5,200,5]
			// lfa: [50]
		},

		N: 1000,
		nyquist: 1,
		steps: 200
	}, function (ran) {
		ran.pipe(process.stdout);
	});
	
### unsupervised and supervised learning mode
	var ran = new RAN({

		learn: function (cb) {
			JSLIB.GET.byStep({
				_Load: [
				]
			}, cb);
		},	

		solve:  {
			batch: 50,
			compress: true,
			interpolate: false,
			lma: [50]
			// bfs: [5,200,5]
			// lfa: [50]
		},

		filter: function (str, ev) {  
			switch (ev.at) {
				case "config":
					str.push(ev);
					break;

				case "batch":
					Log(ev);
					break;

				case "end":
					str.push(ev);
					break;
			}
		},

		K: 2,  // assume 2-state process
		N: 50,  // assume 50 elements in process
	}, function (ran) {
		ran.pipe( function (store) {
			Log(store);
		});
	});
	
## Contributing

See our [issues](https://totem.west.ile.nga.ic.gov/issues.view), [milestones](https://totem.west.ile.nga.ic.gov/milestones.view), [s/w requirements](https://totem.west.ile.nga.ic.gov/swreqts.view),
and [h/w requirements](https://totem.west.ile.nga.ic.gov/hwreqts.view).
	
## License

[MIT](LICENSE)

*/