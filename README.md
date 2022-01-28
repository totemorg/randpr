# RANDPR

Generate or learn various random processes.

In its generation mode, **RANDPR** will generate a first-order, K-state Markov process given 
its first-order, KxK transition probabilites.  **RANDPR** also computes equilibrium probabilities (if they 
exist), mean recurrence times, time to first absorption, absorption probabilites and the statistical 
auto-covariance function.  **RANDPR** can create ergodic, regular, absorptive processes, Weiner 
(stationary in first increments), and wide-sense stationary processes.

In its learning mode, **RANDPR** produces supervised and unsupervised estimates:
MLEs of the underlying transition probabilities, number of coherence intervals (and related SNR),
and the underlying intensity profile.  (MLE for the Weiner process, e.g. first time to exit, have not 
yet been implemented).

Both discrete- and continious-time models are supported in either forward or reverse mode.  

**RANDPR** can be customized with onStep(), onBatch(), onEnd(), and filter() methods to 
sample metrics in both forward and reverse modes.  Both modes can make use of piped 
streams to minimize memory usage.

**RANDPR** supports the following processes.


## Bayes

	K-state process governed by a prescribed conditional independency network:
	
		eqP: [pr, ...] the K equilibrium probs 
		net: [ {var: probs, ...}, ... ] the conditional dependencies 
		
	or expressed as a DAG:
	
		dag: { ... }
	
## Gillespie

	Inhomogenious K-state process with specified transition probabilties:
	
		states: number of states
		
	where its K^2 transition probs are synthesized using the gillespie model.
	
## Markov

	K-state process with specified transition probabilities:
	
		TxPrs: [ [...], ....] the K^2 (K^2-K independent) transition probs 
		
	or:
	
		states: K
		TxPrs: { from: {to: pr, ... } , ... "from, ..." : "to, ..." }	
		
	where from-to transition probs must be specified to conserve prob, i.e. sum_k TxPrs[n][k] = 1.

## Gauss

	Correlated, stateless random process whose parameters are typically derived (see man) 
	for a process with known correlation intervals M = T/Tc or SNR = sqrt{ M / ( 1 + deltaC / M) }.

		values: [ ... ] pc eigen values  [unitless]
		vectors: [ [... ], ...] pc eigen values	[sqrt Hz]
		ref: reference eigenvalue 
		dim: max pc dimension (M = T/Tc )
		mean: mean count in observation interval T

## Wiener

	Stateless process with moving 2nd moment (but stationary in 1st increments) where:
	
		walks: number of walks at each time step (0 disables)
		
	This still needs a quick debug.  May need to define	false Markov parms to init it.
	
## Ornstein

	Stateless Ornstein-Ulenbeck process with:
	
		theta: 0-pi
		a: sigma/sqrt(2 theta)
		
## Mixing

	Gauss mixing process with specified mu,sigma (mean, covar), or specified snr, cone, mixes, oncov, offcov

## Refs

* www.statslab.cam.ac.uk/~rrw1
* www.stat.yale.edu/~pollard
* people.math.gatech.edu/~randall
* www.stat.berkeley.edu/~pitman
* www.math.dartmouth.edu/~pw

## Installation

Clone **RANDPR** from one of its repos:

	git clone https://github.com/totemstan/randpr
	git clone https://sc.appdev.proj.coe/acmesds/randpr
	git clone https://gitlab.west.nga.ic.gov/acmesds/randpr

Dependent modules:

+ **ENUMS** [WWW](https://github.com/totemstan/enums)  [COE](https://sc.appdev.proj.coe/acmesds/enums)  [SBU](https://gitlab.west.nga.ic.gov/acmesds/enums)  

### Manage 

	npm test [ ? || R1 || R2 || ... ]	# unit test
	npm run [ edit || start ]			# Configure environment
	npm run [ prmprep || prmload ]		# Revise PRM

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

where its configuration keys (
[WWW](http://totem.zapto.org/shares/prm/totem/index.html) 
[COE](https://totem.west.ile.nga.ic.gov/shares/prm/totem/index.html) 
[SBU](https://totem.nga.mil/shares/prm/totem/index.html)
)
follow the **ENUMS** deep copy conventions (
[WWW](https://github.com/totemstan/enum) 
[COE](https://sc.appdev.proj.coe/acmesds/enum) 
[SBU](https://gitlab.west.nga.ic.gov/acmesds/enum)
).

The instance can then be piped:

	ran.pipe()
	
to genenerate a process (if configured in the forward/generate mode), or learn process parameters (if 
configured in the reverse/learning mode).


## Program Reference
<details>
<summary>
<i>Open/Close</i>
</summary>
<a name="module_RANDPR"></a>

## RANDPR
Generates various random processes.

**Requires**: <code>module:stream</code>, <code>module:man</code>  
**Example**  
```js
R2.1 - config methods:

	var ran = new RAN({
		p: [.4],
		//markov: [[0.1, 0.9], [0.1, 0.9]]
		//markov: { states: 3, 0: {1: 0.8, 2: 0.1}, 1: {0: 0.1} }
	});
```
**Example**  
```js
R2.3 - config methods:

	var ran = new RAN({
		emP: {
			dims: [3,3],
			weights: [1,1]
		},
		markov: { states: 9, 0: {1: 0.8, 2: 0.1}, 1: {0: 0.1}, "0,1": { "1,0": .4} }
	});
		
```
**Example**  
```js
R2.4 - config methods:

	var ran = new RAN({
		emP: {
			dims: [2,2,2],
			weights: [1,1,1]
		},
		markov: "random"
	});
	
```
**Example**  
```js
R3 sync pipe with various textbook examples, custom filtering with supervised learning validation:
	
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

		batch: 50,  // supervised learning every 50 steps

		filter: function (str, ev) {  
			switch (ev.at) {
				case "config":
					//Trace(ev);
					str.push(ev);
					break;

				case "batch":
					//Trace(ev.s,ev.rel_txpr_error);
					Trace(ev);
					break;

				case "end":
					Trace(ev);
					var
						A = ev.stats.mle_tr_probs,
						B = ev.stats.mle_em_probs,
						H = ev.stats.mle_holding_times;

					Trace("MLEs", {
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
			Trace(store);
		});
		
```
**Example**  
```js
R3.1 - gen process for R3.2 with async pipe to stdout:
	
	var ran = new RAN({

		markov: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3

		batch: 800,  // supervised learning every 50 steps

		N: 1000, 
		filter: function (str,ev) {
			switch (ev.at) {
				case "batch":
				case "config":
				case "end":
					Trace(JSON.stringify(ev));
			}
		},
		steps: 800  
	});

	ran.pipe(process.stdout);
		
```
**Example**  
```js
R3.2 - gen process for R3.3 using async pipe to stdout:

	var ran = new RAN({

		markov: [[0.1, 0.9], [0.1, 0.9]],  // pg142 ex3

		//batch: 50,  // supervised learning every 50 steps

		N: 10,
		//keys: {state:"u", index: "n"},
		filter: function (str,ev) {
			switch (ev.at) {
				case "jump":
					Trace(ev);
					break;
				default:
			}
		},
		steps: 20
	});
	
	ran.pipe(process.stdout);  // stdout evs used in R3.3
		
```
**Example**  
```js
R3.3 - supervised learning with R3.2 evs using sync pipe to store:

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
					Trace( evs ? ` supervising ${evs.length} events` : " supervised" );

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
						Trace(ev);
						str.push(ev);
						break;

					case "batch":
						//Trace(ev.s,ev.rel_txpr_error);
						Trace(ev);
						break;

					case "end":
						//Trace(ev);
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
		Trace(store);
	});
```
<a name="module_RANDPR..filter"></a>

### RANDPR~filter()
Output event filter
				filter: function (str, ev, ran) { // event ev for stream/store str
						switch ( ev.at ) {   // streaming plugins provide an "at" to filter events on
							case "...":
							case "...":
								str.push(ev);	// return the event
						}
					}

**Kind**: inner method of [<code>RANDPR</code>](#module_RANDPR)  
</details>

## Contacting, Contributing, Following

Feel free to 
* submit and status **TOTEM** issues (
[WWW](http://totem.zapto.org/issues.view) 
[COE](https://totem.west.ile.nga.ic.gov/issues.view) 
[SBU](https://totem.nga.mil/issues.view)
)  
* contribute to **TOTEM** notebooks (
[WWW](http://totem.zapto.org/shares/notebooks/) 
[COE](https://totem.west.ile.nga.ic.gov/shares/notebooks/) 
[SBU](https://totem.nga.mil/shares/notebooks/)
)  
* revise **TOTEM** requirements (
[WWW](http://totem.zapto.org/reqts.view) 
[COE](https://totem.west.ile.nga.ic.gov/reqts.view) 
[SBU](https://totem.nga.mil/reqts.view), 
)  
* browse **TOTEM** holdings (
[WWW](http://totem.zapto.org/) 
[COE](https://totem.west.ile.nga.ic.gov/) 
[SBU](https://totem.nga.mil/)
)  
* or follow **TOTEM** milestones (
[WWW](http://totem.zapto.org/milestones.view) 
[COE](https://totem.west.ile.nga.ic.gov/milestones.view) 
[SBU](https://totem.nga.mil/milestones.view)
).

## License

[MIT](LICENSE)

* * *

&copy; 2012 ACMESDS