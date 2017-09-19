/**
@class RANDPR
	[SourceForge](https://sourceforge.net) 
	[github](https://github.com/acmesds/debe.git) 
	[geointapps](https://git.geointapps.org/acmesds/debe)
	[gitlab](https://gitlab.weat.nga.ic.gov/acmesds/debe.git)
	
# RANDPR

In its "forward" (aka "simulation") mode, RANDPR will generate a first-order, K-state Markov process given its first-order, KxK
transition probabilites.  In its "reverse" (aka "realtime") mode, RANDPR recovesr a MLE of these transition probabilities.  Both 
discrete- and continious-time models are supported in either mode.  In its forward mode, RANDPR will compute equilibrium 
probabilities (if they exist), mean recurrence times, time to first absorption, absorption probabilites and the 
statistical auto-covariance function.

As such, RANDPR can create ergodic, regular, absorbtive, processes, as well as Weiner (stationary in first
increments) and wide-sense stationary processes.

In addition to transition probabilites, RANDPR accepts a batch size and nyquist parameters, as well as onStep(), onBatch(), 
onEnd(), and filter() methods to sample metrics in both both modes.  Its forward (soon reverse) of RANDPR makes 
use of piped streams to minimize memory usage.

## Databases

None

## Use

Simply require, configure and start RANDPR:

	var RAN = require("randpr");
	var ran = new RAN({
		key: value, 						// set key
		"key.key": value, 					// indexed set
		"key.key.": value,					// indexed append
		OBJECT: [ function (){}, ... ], 	// add OBJECT prototypes 
		Function: function () {} 			// add chained initializer callback
		:
		:
	}, function () {
		ran.pipe(process.stdout);  // pipe the generated process events to the stdout stream, or ...
		ran.pipe( [], function (events) { // store process events into a list
		});
	});

where its configuration keys follow [ENUM copy()](https://github.com/acmesds/enum) conventions and
are described in its [PRM](/shares/prm/randpr/index.html).

### Ex1

### Ex2

## License

[MIT](LICENSE)

*/