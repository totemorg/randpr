// UNCLASSIFIED 

const
	// nodejs modules
	ENV = process.env,
	OS = require("os"),
	FS = require("fs"),
	STREAM = require("stream"),
	HTTP = require("http"),						//< http interface
	HTTPS = require("https"),					//< https interface	  
	CLUSTER = require("cluster"),				//< for tasking
	VM = require("vm"),							//< for getting
	CP = require("child_process"),				//< spawn OS shell commands
	CRYPTO = require("crypto");					//< crypto for SecureLink

	// 3rd party modules
	IMAP = require('imap'),				// IMAP mail receiver
	MAIL = require('nodemailer'),		// MAIL mail sender
	SMTP = require('nodemailer-smtp-transport'),
	MYSQL = require('mysql'),
	NEO4J = require('neo4j-driver');
	//RSS = require('feed'),				// RSS / ATOM news feeder
	//RSS = require('feed-read'), 		// RSS / ATOM news reader

/*
String.prototype.parseURL = function(opts,base) {
	const 
		url = this+"",
		{username,password,hostname,protocol,pathname,search,port,href} = new URL(url,base);

	//Trace("url", URL(url,base) );

	return Copy( opts || {}, {
		auth: username + ":" + password,
		path: pathname + search,
		protocol: protocol,
		host: hostname,
		port: port,
		href: href
	});
};
	
*/

/**
Provides methods to clock, fetch, enumerate, regulate, stream, fetch and get data.  This module 
documented in accordance with [jsdoc]{@link https://jsdoc.app/}.

### Env Dependencies

	FETCH_PASS = password for local fetching cert
	URL_LEXNEX = URL to lexis-nexis service w user/password credentials
	URL_MYSQL = URL to mysql service w user/password credentials
	URL_NEO4J = URL to neo4j service w user/password credentials
	URL_TXMAIL = URL to smtp email service w user/password credentials
	URL_RXMAIL = URL to imap email service w user/password credentials
	URL_RSSFEED = URL to rss service w user/password credentials
	URL_LEXNEX = URL to lexis-nexis service w user/password credentials

@module ENUMS
@author [ACMESDS](https://totemorg.github.io)

@requires [os](https://nodejs.org/docs/latest/api/)
@requires [cluster](https://nodejs.org/docs/latest/api/)
@requires [fs](https://nodejs.org/docs/latest/api/)
@requires [http](https://nodejs.org/docs/latest/api/)
@requires [https](https://nodejs.org/docs/latest/api/)
@requires [vm](https://nodejs.org/docs/latest/api/)
@requires [cp](https://nodejs.org/docs/latest/api/)
@requires [crypto](https://nodejs.org/docs/latest/api/)
@requires [stream](https://nodejs.org/docs/latest/api/)

@requires [mysql](https://www.npmjs.com/package/mysql)
@requires [neo4j-driver](https://www.npmjs.com/package/neo4j-driver)
@requires [nodemailer](https://www.npmjs.com/package/nodemailer)
@requires [nodemailer-smtp-transport](https://www.npmjs.com/package/nodemailer-smtp-transport)

*/
const
	{ 	
		operators, getList, Regulate,
		Start, Copy, Each, typeOf, Trace, Fetch, Log,
	 	isArray, isObject, isString, isFunction, isEmpty,
		mysqlLogin, neo4jLogin, txmailLogin, rxmailLogin,
		mysqlOpts, neo4jOpts, rxmailOpts, txmailOpts
	} = ENUMS = module.exports = {

	mysqlLogin: ENV.URL_MYSQL ? new URL( ENV.URL_MYSQL ) : null,
	neo4jLogin: ENV.URL_NEO4J ? new URL( ENV.URL_NEO4J ) : null,
	txmailLogin: ENV.URL_TXMAIL ? new URL( ENV.URL_TXMAIL ) : null,
	rxmailLogin: ENV.URL_RXMAIL ? new URL( ENV.URL_RXMAIL ) : null,

/**
*/
	mysqlOpts: {
		// connection options
		connectionLimit: 50,	// max number to create "at once" - whatever that means
		acquireTimeout: 600e3,	// ms timeout during connection acquisition - whatever that means
		connectTimeout: 600e3,	// ms timeout during initial connect to mysql server
		queueLimit: 0,  						// max concections to queue (0=unlimited)
		waitForConnections: true,		// queue when no connections are available and limit reached

		// reserved for ...
		threads: 0, 	// connection threads
		pool: null		// connector
	}, 
/**
*/
	neo4jOpts: { 
		disableLosslessIntegers: true 
	},
/**
*/
	rxmailOpts: {
		secure: true,
		//debug: err => { console.warn(ME+">"+err); } ,
		connTimeout: 10000
	},
/**
*/
	txmailOpts: {
	},

/**
*/
	escapeId: MYSQL.escapeId,
		
/**
*/
	escape: MYSQL.escape,

/**
Start the command prompter for the hosting module with callback to initializer.

@param {String} host Name of hosting module
@param {Array} cbs List of initializer callbacks (ctx defines context to run command)
*/
	Start: (host,cbs) => {
		
		cbs = cbs || {};
		if ( typeof cbs != "function" ) cbs = Copy(cbs, (cmd,ctx) => {
			VM.runInContext( cmd, ctx ); 
		});
		
		const
			VM = require("vm"),
			[x,mod,opt] = process.argv,
			ctx = VM.createContext( Copy( cbs, {} ) ),
			{help} = ctx,
			res = ctx[opt];
			  
		// console.log(process.argv, [opt, host, mod, CLUSTER.isMaster]);
		
		if ( CLUSTER.isMaster && mod && mod.endsWith(host+".js") ) {
			if ( res ) res();
			
			switch (opt) {
				case "?": 
					const opts = ["?","$"].concat(Object.keys(ctx)).join("|");
					Trace( `Usage: node ${host} ${opts}` );
					break;

				case "$":
				default:
					if ( typeof cbs == "function" )
						require("repl").start({
							eval: (cmd, CTX, filename, cb) => {
								//console.log(CTX, filename);

								cmd = cmd.replace(/\n/g,"");

								if ( !cmd ) 
									return;

								else
								if ( cmd.startsWith("<exit") ) 
									process.exit();
								
								else
								if ( cmd.startsWith("<") ) 
									try {
										cbs( FS.readFileSync("./"+cmd.substr(1),"utf-8"), ctx );
									}

									catch (err) {
										console.log("no such file", err);
									}

								else
								if ( cmd.startsWith("??") ) 
									console.log(`
??			this help
?KEY		display contents of variable KEY
?			display all variables
?BOOK:CASE	blog notebook-usecase 
?TOPIC		help on topic
<FILE		embed script FILE
<exit		exit 
`);
									
								else
								if ( cmd.startsWith("?") ) {
									const 
										n = cmd.substr(1);

									if ( book = n.match( /(.*)\:\:(.*)/ ) ) {
										const [nb,uc] = book;
									}

									else
									if ( url = n.match( /(.*)\:(.*)/ ) ) {
									}

									else
									if ( !n )
										console.log(ctx);

									else
									if ( n in ctx )
										console.log( ctx[n] );

									else
									if (help)
										Object.keys(help).forEach( open => {
											if ( n in help[open] ) 
												CP.exec( open.parse$({name: n}) );
										});
								}

								else 
									cbs(cmd,ctx);
							},
							prompt: "$> ", 
							useGlobal: true	// seems to ignore - ctx always "global" but never truely *global*
						});
					break;
			}
		}
	},
	
/**
*/
	Log: console.log,

/**
Trace log message and args. 
*/
	Trace: (msg, ...args) => `enums>>>${msg}`.trace( args ),

/**
Configure enums
@param {Object} opts options
*/
	config: opts => {
		if ( opts ) Copy( opts, ENUMS, "." );
		return ENUMS;
	},
	
/**
Test an object x:

	isString(x), isDate(x), isFunction(x), isArray(x), isObject(x)
	isEmpty(x), isNumber(x), isKeyed(x), isBoolean(x), isBuffer(x)
	isError(x)
*/
			
	typeOf: obj => obj.constructor.name,
	isString: obj => typeOf(obj) == "String",
	isNumber: obj => typeOf(obj)== "Number",
	isArray: obj => typeOf(obj) == "Array",
	isKeyed: obj => Object.keys(obj).length ? true : false,
	isObject: obj => typeOf(obj) == "Object",
	isDate: obj => typeOf(obj) == "Date",
	isFunction: obj => typeOf(obj) == "Function",
	isError: obj => typeOf(obj) == "Error",
	isBoolean: obj => typeOf(obj) == "Boolean",
	isBuffer: obj => typeOf(obj) == "Buffer",
	isEmpty: opts => {
		for ( var key in opts ) return false;
		return true;
	},
	
	/*
	getURL: function (url, cb) {
		Fetch( url, buf => {
			try {
				cb( JSON.parse(buf) );
			}

			catch (err) {
				cb( null );
			}
		});
	},
	*/
		
/**
*/
	getList: function (recs, index, stash) {
		if (index)
			switch (index.constructor.name) {
				case "String":
					const
						opts = {};

					index.split("&").forEach( opt => {
						const
							[lhs,rhs] = opt.split("=");
						opts[lhs] = rhs || "";
					});
					return getList(recs, opts, stash );

				case "Function":
					for ( var i=0; i<recs.length; i++) index(i,recs);
					return recs;

				case "Object":
					const
						res = {},
						where = index["!where"],
						tests = [];

					//Trace(index,stash);

					if ( isEmpty(index) ) 
						if ( stash["!all"] ) 
							Object.index(recs[0] || {}).forEach( key => opts[key] = key );

						else
							return recs;

					if ( where ) {
						recs.forEach( rec => {
							tests.push( VM.runInContext(where,Copy(rec,stash)) );
						});

						delete index["!where"];

						Each( index, (tar,src) => {
							const x = res[tar] = [];
							recs.forEach( (rec,idx) => {
								if ( tests[idx] ) x.push( src ? VM.runInContext(src,Copy(rec,stash)) : rec[tar] );
							});
						});
					}

					else
						try {
							Each( index, (tar,src) => {
								const x = res[tar] = [];
								recs.forEach( rec => {
									x.push( src ? VM.runInContext(src,Copy(rec,stash)) : rec[tar] );
								});
							});
						}
						catch (err) {
						}

					return res;
			}

		else
			return recs;
	},

/**
Copy source hash src to target hash tar.  If the copy is shallow (deep = false), a 
Copy({...}, {}) is equivalent to new Object({...}).  In a deep copy,
(e.g. deep = "."), src keys are treated as keys into the target thusly:

	{	
		A: value,			// sets target[A] = value

		"A.B.C": value, 	// sets target[A][B][C] = value

		"A.B.C.": {			// appends X,Y to target[A][B][C]
			X:value, Y:value, ...
		},	

		OBJECT: [ 			// prototype OBJECT (Array,String,Date,Object) = method X,Y, ...
			function X() {}, 
			function Y() {}, 
		... ]

	} 

@param {Object} src source hash
@param {Object} tar target hash
@param {String} deep copy key 
@return {Object} target hash
*/
	Copy: (src,tar,deep) => {

		for (var key in src) {
			var val = src[key];

			if (deep) 
				switch (key) {
					case Array: 
						val.Extend(Array);
						break;

					case String: 
						val.Extend(String);
						break;

					case Date: 
						val.Extend(Date);
						break;

					case Object: 	
						val.Extend(Object);
						break;

					/*case "Function": 
						this.callStack.push( val ); 
						break; */
						
					default:

						var 
							keys = key.split(deep), 
							Tar = tar,
							idx = keys[0],
							N = keys.length-1;

						//Trace("##########", deep, keys, N);
						
						for ( var n=0; n < N ;  idx = keys[++n]	) { // index to the element to set/append
							if ( idx in Tar ) {
								if ( !Tar[idx] ) Tar[idx] = new Object();
								Tar = Tar[idx];
							}

							else
								Tar = Tar[idx] = new Object(); //new Array();
						}

						if (idx)   	// not null so update target
							Tar[idx] = val;
						
						else  		// null so append to target
						if (val.constructor == Object) 
							for (var n in val) 
								Tar[n] = val[n];
						
						else
							Tar.push( val );
				}

			else
				tar[key] = val;
		}

		return tar;
	},

/**	 
Enumerate Object A over its keys with callback cb(key,val).

@param {Object} A source object
@param {Function} cb callback (key,val) 
*/
	Each: ( A, cb ) => {
		Object.keys(A).forEach( key => cb( key, A[key] ) );
	},
	
/*
Stream a src array, object or file using:

	Stream(src, opts, (rec,key,res) => {
		if ( res ) // still streaming 
			res( msg || undefined )  // pass undefined to bypass msg stacking

		else 
			// streaming done so key contains msg stack
	})

@param {Object | Array | String} src source object or array
@param {Function} cb( rec || null, key, res ) Callback	

@example
Serialize a list:

	function fetcher( rec, info => { 
	});

	[ rec, ...].serial( fetcher, (rec, fails) => {
		if ( rec ) 
			// rec = record being serialized
		else
			// done. fails = number of failed fetches
	}

@example
### Serialize a string:

	function fetcher( rec, ex => {
		// regexp arguments rec.arg0, rec.arg1, rec.arg2, ...
		// rec.ID = record number being processed
		return "replaced string";
	});

	"string to search".serial( fetcher, regex, "placeholder key", str => { 
		// str = final string with all replacements made
	});

*/
	/*
	Stream: (src,opts,cb) => {
		
		var
			msgs = []; 
		
		if ( isString(src) ) 	// stream file
			src.streamFile( opts, recs => {
				if ( recs )
					recs.forEach( (rec,idx) => {
						cb( rec, idx, msg => {
							if ( msg != undefined ) 
								msgs.push( msg );
						});
					});
				
				else 	// signal end
					cb(null, msgs);
			});
		
		else {	// stream array/object
			var 
				returns = 0, 
				calls = src.forEach ? src.length : Object.keys(src).length;

			Each( src, (key, rec) => {
				cb( rec, key, msg => {
					if ( msg != undefined ) 
						msgs.push( msg );
						/ *
						if ( A.forEach ) 
							msgs.push( msg );
						else
							msgs[key] = msg;
						* /

					if ( ++returns == calls ) // signal end
						cb( null, msgs );
				});
			});

			if ( !calls ) cb( null, msgs );
		}
	}, */

	Clock: Clock,

	// Fetch options

	//defHost: ENV.SERVICE_MASTER_URL || "http://UNDEFINED:UNDEFINED",
	
/**
Fetch quick SITEREFs

	https://www.programmableweb.com/search/military

	ACLED
	https://www.programmableweb.com/api/acled-rest-api-v25
	ACCT teliy40602@plexfirm.com / ACLEDsnivel1
	API https://developer.acleddata.com/rehd/cms/views/acled_api/documents/API-User-Guide.pdf
	SITE https://developer.acleddata.com/
	The Armed Conflict Location & Event Data Project (ACLED) is a real-time data and and crisis analysis and mapping project on political violence and protest across the world. ACLED's mission is to produce dis-aggregated, locally informed data and analysis globally in real time. An ACLED REST API enables users to retrieve data about actors, actor type, country, region and get data in JSON, XML, CSV or text. Filter data by date, area, pagination, etc.

	Animetrics FIMS
	https://www.programmableweb.com/api/animetrics-fims-cloud-rest-api
	http://animetrics.com/fims-cloud
	Aimed at the law enforcement, security and military intelligence industries, Animetrics' FaceR Identity Management Solution (FIMS) allows organizations to bring mobile security and video surveillance facial-biometrics applications into the field for use in real time. FIMS Cloud is used to centralize and access a user's cloud based photographic stores using facial recognition. FIMS Cloud organizes, searches and centralizes access to photographic stores using 1:many web service based verification engine. Access to the service is provided via a RESTful API. Public documentation is not available.

	Navlost WXT Weather Tesseract
	https://www.programmableweb.com/api/navlost-wxt-weather-tesseract
	The WXT Weather Service provides atmospheric weather information through a REST architecture, HTTP requests, and JSON formats. It integrates METAR/TAF information, sun, and moon calculations, targeting aviation and energy applications. Advanced features include: -Upper atmosphere information (e.g., research, aviation, rocketry, military) -Automated, push-type notification of arbitrary weather-related events (alert service) -Calculation of arbitrary results derived from weather forecast information via a server-side scripting language. The default response type is application/json, although other formats may be supported. At the present time, there is partial support for comma-separated value (CSV) responses.
	https://wxt.navlost.eu/api/v1/
	https://wxt.navlost.eu/doc/api/
*/
	sites: {
		"*acled": "https://api.acleddata.com/{data}/{command}",
		"*navlost": "https://wxt.navlost.eu/api/v1/${find}",
		"*nmi": "http://api.met.no/weatherapi/?opt=${x}"
	},
/**
Max files to Fetch when indexing a folder
*/
	maxFiles: 1000,
/**
Fetch wget/curl maxRetry	
*/
	maxRetry: 5,		
/**
Legacy Fetching certs
*/
	certs: { 			
		//pfx: FS.readFileSync(`./certs/fetch.pfx`),
		//key: FS.readFileSync(`./certs/fetch.key`),
		//crt: FS.readFileSync(`./certs/fetch.crt`),
		//ca: "", //FS.readFileSync(`./certs/fetch.ca`),			
		_pfx: `./certs/fetch.pfx`,
		_crt: `./certs/fetch.crt`,
		_key: `./certs/fetch.key`,
		_ca: `./certs/fetch.ca`,
		_pass: ENV.FETCH_PASS
	},
	
/**
Regulate a task defined by options `opts`

	every 	= "NUM sec||min||hr||... "
	start	= DATE  
	end		= DATE  
	on		= NUM  
	off		= NUM  
	util	= NUM  
	batch	= INT  
	limit	= INT  

with callbacks to

	taskcb( recs ) to process a `recs`-batch
	feedcb( step ) to feed the queue via `step(recs)` 

When a `feedcb` is provided, the mandatory `taskcb` is placed into a 
buffered workflow that terminates when the `recs`-batch goes null.  

If no `feedcb` is provided, the `taskcb` is periodically executed (i.e.
in an unbuffered workflow) using null `recs`-batches.

The regulated task is monitored/managed by the supplied options (defaults)

	task 	= notebook being regulated ("task")
	name	= usecase being regulated ("case")
	client	= task owner ("system")
	watch	= QoS task watchdog timer [s]; 0 disables (60)

A nonzero QoS sets a tasking watchdog timer to manage the task.  A credit
deficient client is signalled by calling `feedcb(null)`.

To establish the task as a proposal, set Sign0 = 1 in the taskDB: in so 
doing, if Sign1 , ... are not signed-off (eg not approved by a task oversight
commitee) before the proposal's start time, the task will be killed.

Regulate uses the following DBs:

	userDB = openv.profiles client credit/billing information
	taskDB = openv.queues tasking/billing information
	snapDB = openv.<task> holds the task snapshot context 

@param {Object} opts Task regulation options hash
@param {Function} taskcb(recs) Process record batch recs 
@param {Function} [feedcb(step)] Feed a record batch using step(recs)
@returns {Clock} Clock built for regulation options
*/
	Regulate: (opts,taskcb,feedcb) => {
		// ==================== configuration options

		const
			snapStore = "", 		// "Save_snap" to enable snapshots
			userDB = "openv.profiles",
			taskDB = "openv.queues",
			snapDB = "app.??";

		// ===================== configuration options

		function cpuUtil() {				// return average cpu utilization
			var avgUtil = 0;
			var cpus = OS.cpus();

			cpus.forEach( cpu => {
				idle = cpu.times.idle;
				busy = cpu.times.nice + cpu.times.sys + cpu.times.irq + cpu.times.user;
				avgUtil += busy / (busy + idle);
			});
			return avgUtil / cpus.length;
		}

		function memUtil() {
			return OS.freemem() / OS.totalmem();
		}

		function manageTask( sql, task ) {
			if ( task )
				if ( task.Kill ) 
					stopTask(sql, task.ID, true);
				else {
				}

			else // someone deleted the task
				stopTask( sql, 0, true );
		}

		function startTask( sql, cb ) {
			sql.query(
				"SELECT Credit FROM ?? WHERE Client=?", [userDB, Client], (err,profs) => {

				if ( err ) 
					Trace("cant start task - check if mysql db started");
					
				else
				if ( prof = profs[0] )
					if ( prof.Credit )
						sql.query(  			// add job to the job queue
							"INSERT INTO ?? SET ? ON DUPLICATE KEY UPDATE Run=Run+1,?",
							[taskDB,queueRec,queueRec],
							(err,info) => {  	// increment work backlog for this job

							const
								taskID = err ? 0 : info.insertId;

							if ( taskID && watchdog ) 	// setup task watchdog
								watchTask = setInterval( jobs => {
									sqlThread( sql => {
										sql.query("SELECT * FROM ?? WHERE ID=?", [taskDB,taskID], (err,recs) => manageTask( sql, recs[0] ) );
									});
								}, watchdog, jobs );

							cb( taskID );
						});

					else	// signal client has no credit 
						cb( 0 );

				else	// signal no such client
					cb( 0 );
			});
		}

		function stopTask( sql, taskID, kill ) {
			if ( kill ) {								// kill all pending jobs and flush queue
				jobs.forEach( (job,idx) => {			// clear jobs
					clearTimeout( job.timer );
					delete jobs[ idx ];
				});
				jobs.length = 0;						// flush queue
			}

			if (watchTask) clearTimeout(watchTask);		// stop watching the task

			sql.query( 
				"UPDATE ?? SET Departed=now(),Finished=1 WHERE ID=?", 
				[taskDB,taskID] );

			sql.query(  			// update clients credit
				"UPDATE ?? SET Credit=Credit-? WHERE Client=?",
				[ userDB, jobsTasked, Client ] );
		}

		function updateTask( sql, taskID ) {
//Trace("update", taskID, evsTasked, jobsTasked);
			sql.query(  			// update job queue
				"UPDATE ?? SET " +
				"Done=Done+1, " +
				"Age=datediff(now(),Arrived)*24, " + 
				"State=Age/Done, " +
				"ECD=?,Snaps=Snaps+1,Events=?,Batches=? WHERE ID=?", 
				[ taskDB, clock.next, evsTasked, jobsTasked, taskID ] );		// "ECD=date_add(Arrived, interval State*Work hour), "
		}
		
		function traceTask( msg ) {
			`${Task}-${Name}`.trace( msg );
		}

		const
			{ round } = Math,
			{ every,start,end,on,off,util,watch,limit,batch,name,client,task } = opts,
			{ Every,Start,End,On,Off,Util,Watch,Limit,Batch,Name,Client,Task } = Opts = {
				// event clocking
				
				Every: every+"",					// N[sec||min||...]
				Start: start ? new Date( start ) : new Date(),
				End: end || null,
				On: round(on || off * util || 0),						// [steps]
				Off: round( off || (util ? on * (1-util)/util : 0) ),	// [steps]
				//Cycle: on + off,										// [steps]
				Util: util || on / (on+off),							// [0:1]

				// file regulation
				
				Batch: batch || 0,										// [recs]
				Limit: limit || 0,										// [recs]

				// job regulation
				
				Watch: watch || 60,										// [secs]
				Task: task || taskcb.name || "task",
				Name: name || "case",
				Client: client || "system"
			},
			clock = new Clock( "",Every,On,Off,Start,End );

//Trace( "Regulating", clock);
		
		const 
			jobs = [],
			queueRec = {
				// mysql unique keys should not be null
				Client: Client,
				Class: "fetch",
				Task: Task,
				Name: Name,

				// job qos
				QoS: Watch,
				Priority: 0,

				// others 
				Notes: [
					"Notebook".link( `/${Name}.run` ), 
					"Client: "+Client.link( `mailto:${client}?subject=hello` ),
					// ((profile.Credit>0) ? "funded" : "unfunded").tag( url ),
					"RTP".link( `/rtpsqd.view?task=${Name}` ),
					"PMR".link( `/briefs.view?options=${Name}` ),
					//"Source".link(pipePath),
					"Metrics".link( `/plot.view?w=800&h=600&x=Save_baseline$t&y=Save_baseline$y&src=/${Name}` )
				].join(" || "),

				Arrived	: new Date(),
				Departed: null,
				Classif : "(U)",

				// flags
				Kill: 0,
				Billed: 0,
				Flagged: 0,
				Finished: 0,
				Funded: 1,

				// signoffs
				Sign0: 0,
				Sign1: 0,
				Sign2: 0,
				Sign3: 0,

				// metrics
				Events: 0,
				Batches: 0,
				Snaps: 0,

				// Completion estimates
				Age: 0,
				State: 0,
				ECD: null,
				cpuUtil: cpuUtil(),
				memUtil: memUtil(),
				Work: 0,
				Done: 1
			},
			watchdog = Watch*1e3 || 0;	// task watchdog interval in msecs (0 disables)

		var							// task metrics
			watchTask = 0,			// watchdog timer
			evsTasked = 0,			// # of records routed to taskcb
			jobsTasked = 0,			// # of batches routed to taskcb
			evsQueued = 0,			// # of records provided by feedcb
			jobsQueued = 0;			// # of batches queued by feedcb

		sqlThread( sql => {
			startTask( sql, taskID => {
				function stepTask( idx ) {

					function startSnap( cb ) {	// cb(sql,ctx)
						sqlThread( sql => {
							sql.query(		// check approval status
								"SELECT Sign0,Sign1,Sign2,Sign3 FROM ?? WHERE ID=?",
								[ taskDB, taskID ],
								(err,tasks) => {

//Trace("start snap", err, tasks);
								if ( task0 = tasks[0] ) {
									const {Sign0,Sign1,Sign2,Sign3} = task0;

									if ( Sign0 && (!Sign1 || !Sign2 || !Sign3) ) {	// signers have rejected task proposal so ...
										stopTask( sql, taskID, true );		// kill task
										cb( sql, null );					// signal task dead
									}

									else
									if ( snapStore ) {	// provide snapshot context
										//sql.query("START TRANSACTION");
										sql.query(
											`SELECT ${snapStore} AS ctx FROM ?? WHERE Name=? LIMIT 1`, // + "FOR UPDATE", 
											[snapDB,Name], 
											(err,snaps) => {

//Trace("load snapshot", err, snaps);
												if ( err ) 	// signal nonexitant task
													cb( sql, null );

												else
												if ( snap = snaps[0] )
													try {
														cb( sql, JSON.parse( snap.ctx ) );
													}
												
													catch (err) {
														cb( sql, null );
													}

												else	// signal usecase does not exist
													cb( sql, null );
											});
									}
									
									else	// provide dummy context
										cb(sql, {});
								}

								else  	// signal task deleted
									cb( sql, null );
							});
						});
					}

					function saveSnap( sql, save ) {
						const 
							savectx = {
								Host: Task,
								Name: Name
							};

						//Trace(".............save snap", savectx);
						savectx[snapStore] = save;
						Each( save, (key,val) => {
							//Trace(".........dump?", key);
							if ( key.startsWith("_") ) {
								//Trace(".............special dump key", key);
								savectx[key] = val;
								delete save[key];
							}
						});
						sql.saveContext(savectx);
					}

					function stopSnap( sql, ctx, kill ) {
						stopTask( sql, taskID, kill );
						/*taskcb( null, ctx, save => {
							if ( save ) 
								saveSnap(sql,save);
						});*/
					}

//Trace("step task", idx);

					startSnap( (sql,ctx) => {	// load previous snapshot context

//Trace("====step sql", idx, sql?true:false);

						if ( ctx ) {			// valid context so step the task
							const
								recs = jobs[idx].recs;

							if ( recs ) {	// advance metrics 
//Trace("====step metrics", idx, "recs", recs.length);

								evsTasked += recs.forEach ? recs.length : 1;
								//jobsTasked ++;
								// traceTask( `drained ${jobsTasked} batches ${evsTasked} events` );
							}

							traceTask( `completed ${Math.trunc(100*jobsTasked/jobs.length)}%` );
							
							taskcb( recs );
							
							if ( ++jobsTasked >= jobs.length ) 		// no more jobs so issue soft halt
								stopTask( sql, taskID, false );
							
							else
								updateTask( sql, taskID );
							
							delete jobs[idx];		// free job records

							/*
							if ( save ) 	// update and reset snapshot metrics
								saveSnap(sql,save);

							if ( clock.end && (new Date() >= clock.end) )  		// job past its PoP so issue hard halt
								stopSnap( sql, ctx, true );

							if ( ++jobsTasked >= jobs.length ) 		// no more jobs so issue soft halt
								stopSnap( sql, ctx, false);

							else*/

							//sql.query("COMMIT");
						}

						else 					// no context (task rejected) so kill the task
							traceTask( "deleted" );
							//stopTask( sql, taskID, true );

					});
				}

				if ( taskID )		// task succesfully started
					if ( feedcb ) 	// requesting buffered task
						feedcb( recs => {	// get a batch of records
							if ( recs ) {		// queue this record batch
								const
									Recs = recs.forEach ? [] : recs;	// retain either the recs batch or the recs file path

//Trace("feeding", recs.forEach?recs.length : recs);

								if ( recs.forEach ) 		// clone record batch
									recs.forEach( rec => Recs.push( new Object(rec) ));

//Trace("PUSH STEP", jobs.length);
								jobs.push({				// add a job to the task queue
									recs: Recs,			// hold batch at this step
									idx: jobs.length,	// at step idx=job.length=#feeds
									timer: 	setTimeout( 
										idx => {
//Trace("RUN STEP", idx); 
											stepTask( idx );
										}, 

										clock.tick( wait => {	// add a snapshot signal 

//Trace("PUSH SNAP", jobs.length, wait);
											jobs.push({	// callback notebook with null batch to signal snapshot
												recs: null,
												idx: jobs.length,
												timer: setTimeout( idx => {
//Trace("RUN SNAP", idx);
														stepTask( idx );
													}, 

													wait, 

													jobs.length )
											});	

										}), 

										jobs.length )
								});

								jobsQueued++;
								evsQueued += recs.forEach ? recs.length : 1;
								// traceTask( `queued ${jobsQueued} batches ${evsQueued} events` );
							}

							else {	// prime task and schedule end-of-stream job
								jobs.push({				// add a job to the task queue
									recs: null,
									idx: jobs.length,
									timer: setTimeout( 
										idx => {
//Trace("BATCH STEP", idx); 
											stepTask( idx );
										}, 

										clock.tick( ), 

										jobs.length )
								});
								sql.query(
									"UPDATE ?? SET Work=?,ECD=date_add(Arrived, interval State*Work hour) WHERE ID=?", 
									[taskDB, jobs.length,taskID] );
							}
						} );

					else {			// requesting unbuffered task
						setTimeout( 
							idx => {		// initial delay to syncup on requested interval
								jobs.push({
									recs: null,
									idx: jobs.length,
									timer: setInterval( 
											idx => {
//Trace("unreg task", idx);
												taskcb( null );
												updateTask(sql,taskID);
/*												, jobs[idx], save => {
//Trace("unreg job", clock.tick(), clock.every);
													updateTask(sql,taskID);
												}); */
											},  

											clock.tick(), 

											jobs.length )
								});
							}, 

							clock.tick( ), 

							0);
					}

				else				// task could not be started
				if ( feedcb ) 		// signal error condition
					feedcb( null );

				else				// log error
					Log( new Error(`client ${Client} cant start`) );
			});
		} );
		
		return clock;
	},
		
/**
GET (PUT || POST || DELETE) information from (to) a `ref` url

	PROTOCOL://HOST/FILE ? QUERY & FLAGS
	SITEREF

given a `cb` callback function (or a `data` Array || Object || null).

The `ref` url specifies a PROTOCOL

	http(s) 	=	http (https) protocol
	curl(s) 	=	curl (curls uses certs/fetch.pfx to authenticate)
	wget(s)		=	wget (wgets uses certs/fetch.pfx to authenticate)
	mask 		=	http access via rotated proxies
	file		=	file or folder path
	notebook	=	selected notebook record
	lexnex 		=	Lexis-Nexis oauth access to documents

When a file path is "/"-terminated, a folder index is returned.  File paths 
may also contain wild-* cards.  

Use the FLAGS

	_every 	= "NUM sec||min||hr||..."
	_start	= DATE  
	_end	= DATE  
	_on		= NUM  
	_off	= NUM  						
	_util	= NUM  

to regulate the fetch in a job queue with callbacks to the `cb` task.  Use
the FLAGS

	_watch	= NUM  
	_task 	= "job task name"
	_name	= "job case name"
	_client = "job owner"

to monitor the task in the job queues.  Use the FLAGS

	_batch	= NUM
	_limit	= NUM
	_keys	= [...]
	_comma	= "delim"
	_newline= "delim"

to read a csv-file and feed record batches to the `cb` callback. 

@param {String} ref source URL
@param {string | array | function | null} [cb] callback or data 
@param {function} [cb] optional callback when first cb is data

@example
Fetch( ref, text => {			// get request
})

@example
Fetch( ref, [ ... ], stat => { 	// post request with data hash list
})

@example
Fetch( ref, { ... }, stat => { 	// put request with data hash
})

@example
Fetch( ref, null, stat => {		// delete request 
})

*/
	Fetch: (ref, data, cb) => {	//< data fetcher

		function sha256(s) { // reserved for other functionality
			return CRYPTO.createHash('sha256').update(s).digest('base64');
		}

		function request(proto, opts, data, cb) {
//Trace(">>>fetch req opts", opts);
			//delete opts.auth;
			
			const Req = proto.request(opts, Res => { // get reponse body text
				var body = "";
				Res.on("data", chunk => body += chunk.toString() );

				Res.on("end", () => {
					//Trace('fetch statusCode:', Res.statusCode);
					//Trace('fetch headers:', Res.headers['public-key-pins']);	// Print the HPKP values

					if ( cb ) 
						cb( body );

					else
						data( body );
				});
			});

			Req.on('error', err => {
				Trace("fetch error", err);
				(cb||data)("");
			});

			switch (opts.method) {
				case "DELETE":
				case "GET": 
					break;

				case "POST":
				case "PUT":
					//Trace("fetch post", data);
					Req.write( data ); //JSON.stringify(data) );  // post parms
					break;
			}					

			Req.end();
		}

		function agentRequest(proto, opts, sql, id, cb) {
			var 
				body = "",
				req = proto.get( opts, res => {
					var sink = new STREAM.Writable({
						objectMode: true,
						write: (buf,en,sinkcb) => {
							body += buf;
							sinkcb(null);  // signal no errors
						}
					});

					sink
					.on("finish", () => {
						var stat = "s"+Math.trunc(res.statusCode/100)+"xx";
						Trace("agent body", body.length, ">>stat",res.statusCode,">>>",stat);

						sql.query("UPDATE openv.proxies SET hits=hits+1, ?? = ?? + 1 WHERE ?", [stat,stat,id] );

						cb( (stat = "s2xx") ? body : "" );
					})
					.on("error", err => {
						Trace("agent error", err);
						cb("");
					});

					res.pipe(sink);
				});

			req.on("socket", sock => {
				sock.setTimeout(2e3, () => {
					req.abort();
					Trace("agent timeout");
					sql.query("UPDATE openv.proxies SET hits=hits+1, sTimeout = sTimeout+1 WHERE ?", id);
				});

				sock.on("error", err => {
					req.abort();
					Trace("agent refused");
					sql.query("UPDATE openv.proxies SET hits=hits+1, sRefused = sRefused+1 WHERE ?", id);
				});
			});

			req.on("error", err => {
				Trace("agent abort",err);
				sql.query("UPDATE openv.proxies SET hits=hits+1, sAbort = sAbort+1 WHERE ?", id);
			});
		}
		
		function oauthRequest(opts, {grant, token ,doc}, res) {
			
			token.method = "POST";
			token.headers = {
				"Content-Type": "application/x-www-form-urlencoded"
			};
			
			request(HTTPS, token, grant, token => {		// request access token
//Trace("token", token);
				try {
					const 
						Token = JSON.parse(token);

					opts.protocol = "https:";
					opts.headers = {
						Authorization: Token.token_type + " " + Token.access_token,
						Accept: "application/json;odata.metadata=minimal",
						Host: "services-api.lexisnexis.com",
						Connection: "Keep-Alive",
						"Content-Type": "application/json"
					};
					//delete opts.auth;

Trace("oauth token request", opts );
					request(HTTPS, opts, search => {	// request a document search
						if ( docopts = doc ) 			// has document fetch url so ...
							try {	
								const
									Search = JSON.parse(search),
									rec = Search.value[0] || {},
									doclink = rec['DocumentContent@odata.mediaReadLink'];

								//Trace( Object.keys(Search) );
								if ( doclink ) {
									if (false)
										Trace("get doc", {
											doclink: doclink , 
											href: doc.href, 
											reckeys: Object.keys(rec), 
											ov: rec.Overview, 
											d: rec.Date
										});

									if ( docopts.length ) // string so do fetch
										Fetch( docopts + doclink, doc => {
											res(doc);
										});

									else {
										docopts.pathname += doclink;		//< check enums_lexnex code
										docopts.headers = {
											"Content-Type": "application/x-www-form-urlencoded"
										};
										docopts.method = "GET";
										Trace("doc request", docopts);
										request( HTTPS, docopts, doc => res(doc) );
									}
								}
							}

							catch (err) {
								Trace("oauth bad search",err);
							}

						else							// just doing a search so ...
							res( search );
					});
				}

				catch (err) {
Trace("oauth bad token", token);
					res(null);
				}
			});
		}
		
		function fetch(res) {	// res([...] || {...} || "...")
			
			function matchFiles( path, cb ) {
				const
					[ x, folder, file ] = path.match( /(.*)\/(.*)/ ),
					search = new RegExp( file.replace( /\*/g, ".*" ) );

				if ( file.indexOf("*") >= 0 )
					Fetch( "file:"+folder+"/", files => {
						files.forEach( name => {
							if ( name.match(search) ) 
								cb( folder+"/"+name );
						});
					});
				
				else
					cb( path );
			}
			
			switch ( opts.protocol ) {
				case "curl:": 
					CP.exec( `curl --retry ${maxRetry} ` + path.replace(opts.protocol, "http:"), (err,out) => {
						res( err ? "" : out );
					});
					break;

				case "curls:":
					CP.exec( `curl --retry ${maxRetry} -gk --cert ${cert._crt}:${cert._pass} --key ${cert._key} --cacert ${cert._ca}` + url.replace(protocol, "https:"), (err,out) => {
						res( err ? "" : out );
					});	
					break;

				case "wget:":
					CP.exec( `wget --tries=${maxRetry} -O ${wout} ` + path.replace(opts.protocol, "http:"), err => {
						res( err ? "" : "ok" );
					});
					break;

				case "wgets:":
					CP.exec( `wget --tries=${maxRetry} -O ${wout} --no-check-certificate --certificate ${cert._crt} --private-key ${cert._key} ` + wurl.replace(protocol, "https:"), err => {
						res( err ? "" : "ok" );
					});
					break;

				case "https:":
					/*
					// experiment pinning tests
					opts.checkServerIdentity = function(host, cert) {
						// Make sure the certificate is issued to the host we are connected to
						const err = TLS.checkServerIdentity(host, cert);
						if (err) {
							Trace("tls error", err);
							return err;
						}

						// Pin the public key, similar to HPKP pin-sha25 pinning
						const pubkey256 = 'pL1+qb9HTMRZJmuC/bB/ZI9d302BYrrqiVuRyW+DGrU=';
						if (sha256(cert.pubkey) !== pubkey256) {
							const msg = 'Certificate verification error: ' + `The public key of '${cert.subject.CN}' ` + 'does not match our pinned fingerprint';
							return new Error(msg);
						}

						// Pin the exact certificate, rather then the pub key
						const cert256 = '25:FE:39:32:D9:63:8C:8A:FC:A1:9A:29:87:' + 'D8:3E:4C:1D:98:JSDB:71:E4:1A:48:03:98:EA:22:6A:BD:8B:93:16';
						if (cert.fingerprint256 !== cert256) {
							const msg = 'Certificate verification error: ' +
							`The certificate of '${cert.subject.CN}' ` +
							'does not match our pinned fingerprint';
							return new Error(msg);
						}

						// This loop is informational only.
						// Print the certificate and public key fingerprints of all certs in the
						// chain. Its common to pin the public key of the issuer on the public
						// internet, while pinning the public key of the service in sensitive
						// environments.
						do {
							console.log('Subject Common Name:', cert.subject.CN);
							console.log('  Certificate SHA256 fingerprint:', cert.fingerprint256);

							hash = crypto.createHash('sha256');
							console.log('  Public key ping-sha256:', sha256(cert.pubkey));

							lastprint256 = cert.fingerprint256;
							cert = cert.issuerCertificate;
						} while (cert.fingerprint256 !== lastprint256);

						};
					*/
					/*
					opts.agent = new HTTPS.Agent( false 
						? {
								//pfx: FS.readFileSync(cert._pfx),	// pfx or use cert-and-key
								cert: FS.readFileSync(cert._crt),
								key: FS.readFileSync(cert._key),
								passphrase: cert._pass
							} 
						: {
							} );
						*/
					opts.rejectUnauthorized = false;
					request(HTTPS, opts, data, res);
					break;

				case "http:":
//Trace("http", opts);

					request(HTTP, opts, data, res);
					break;

				case "mask:":
				case "mttp:":	// request via rotating proxies
					opts.protocol = "http:";
					sqlThread( sql => {
						sql.query(
							"SELECT ID,ip,port FROM openv.proxies WHERE ? ORDER BY rand() LIMIT 1",
							[{proto: "no"}], (err,recs) => {

							if ( rec = recs[0] ) {
								opts.agent = new AGENT( `http://${rec.ip}:${rec.port}` );
								Trace("agent",rec);
								agentRequest(HTTP, opts, sql, {ID:rec.ID}, res );
							}
						});
					});
					break;

				case "book:":
				case "notebook:":
					sqlThread( sql => {
						sql.query( Type
							? "SELECT * FROM app.? WHERE Name=?"
							: "SELECT * FROM app.?", 

							[ Name, Type ], 

							(err,recs) => res( err ? "" : JSON.stringify(recs) ) );
					});
					break;

				case "file:":	// requesting file or folder index
//Trace("index file", [path], opts);
					matchFiles( "."+Path, src => {	
						const recs = [];

						switch (Type) {
							case "/":
								try {
									FS.readdirSync(src).forEach( file => {
										var
											ignore = file.startsWith(".") || file.startsWith("~") || file.startsWith("_") || file.startsWith(".");

										if ( !ignore && recs.length < maxFiles ) 
											recs.push( (file.indexOf(".")>=0) ? file : file+"/" );
									});
									res( recs );
								}

								catch (err) {
									res( null );
								}
								break;

							case "json":
								try {
									FS.readFile(src, "utf8", (err,buf) => res( err ? null : JSON.parse(buf) ) );
								}

								catch (err) {
									res( null );
								}
								break;

							case "null":
								const 
									dummy = {},
									{batch,limit} = flags;
								var
									read = 0;

								for ( var n=0; n<batch; n++ ) recs.push( dummy );

								while ( read < limit ) { 
									cb(recs);
									read += batch;
								}

								res(null);
								break;

							case "csv":	
								flags.keys = [];
//Trace("stream csv", flags,query);

								src.streamFile( flags, recs => {
	//Trace("csv batch", recs ? recs.length : null);
									res( recs ? recs.get(query) : null );
								});
								break;

							case "list":
								flags.keys = buf => buf;
								Copy(query,flags);

								src.streamFile( flags, paths => {
									paths.forEach( path => res( path ) );
								});
								break;

							case "stream":
							case "str":
								src.streamFile( flags, recs => res( recs ? recs.get(query) : null ) );
								break;

							case "jpg":							
							case "png":							
							case "ico":							
								FS.readFile(src, "buffer", (err,buf) => res( err ? null : buf ) );
								break;
							
							case "txt":
							default:
								FS.readFile(src, "utf8", (err,buf) => res( err ? null : buf ) );
								break;
						}
					});
					break;

				case "lex:": 	// lexis-nexis search only
				case "lexnexsearch":
					oauthRequest(opts, {
						grant: "grant_type=client_credentials",  // &scope=http://auth.lexisnexis.com/all
						token: new URL( URL_LEXNEX_TOKEN )
					}, res );
					break;
					
				case "lexis:": 	// lexis-nexis search and get doc
				case "lexnex":
					oauthRequest(opts, {
						grant: "grant_type=client_credentials",  // &scope=http://auth.lexisnexis.com/all
						token: new URL( URL_LEXNEX_TOKEN ),
						doc: new URL( URL_LEXNEX_GETDOC )
					}, res );
					break;

				default:	
Trace( "fetch bad protocol", opts );
					res(null);
			}
		}
		
		const
			{certs,maxRetry,maxFiles,sites} = ENUMS,
			bufferFiles = true,
			query = {},
			flags = {},
			[Path,Name,Type,Area,opts] = (sites[ref] || ref).parsePath(query,{},flags,{}),
			crud = {
				"Function": "GET",
				"Array": "PUT",
				"Object": "POST",
				"Null": "DELETE"
			},
			cert = certs.fetch,		// for wget-curl
			[wurl,wOut] = opts.pathname.split("////"),
			wopt = wOut || "./temps/wget.jpg",
			res = cb || data || console.log,	// response callback
			method = opts.method = crud[ data ? typeOf(data) : "Null" ] ;

		/*
		// opts.cipher = " ... "
		// opts.headers = { ... }
		// opts.Cookie = ["x=y", ...]
		// opts.port = opts.port ||  (protocol.endsWith("s:") ? 443 : 80);
		if (opts.soap) {
			opts.headers = {
				"Content-Type": "application/soap+xml; charset=utf8",
				"Content-Length": opts.soap.length
			};
			opts.method = "POST";
			opts.pathname = 
				ref => ref.parse$(query) || ref.tag("?",query)
		}*/

Trace("FETCH",ref, query, flags);

		if ( flags.every )	// regulated fetch
			switch (opts.protocol) {
				case "file:":	// regulated fetch
					if (bufferFiles)		// buffered
						Regulate(flags, res, step => {	// feed records to queue
							if (step) 
								fetch(step);

							else
								Trace("client exhausted credit");
						} );
					
					else					// unbuffered
						Regulate(flags, recs => fetch( res ));

					break;
					
				default:		// unbuffered
					Regulate(flags, recs => fetch( res ));
			}

		else				// unregulated fetch
			fetch(res);
	}

};

/**
@module ENUMS.Array
*/
Copy({
/**
*/
	stream: function (cb) {
		const
			calls = this.length,
			msgs = [];
		
		var 
			returns = 0;

		this.forEach( (rec,key) => {
			cb( rec, key, msg => {
				if ( msg != undefined ) 
					msgs[key] = msg;

				if ( ++returns == calls ) // signal end
					cb( null, msgs );
			});
		});

		if ( !calls ) cb( null, msgs );		
	},
	
/**
*/
	Extend: function (con) {
		this.forEach( function (proto) {
			//console.log("ext", proto.name, con);
			con.prototype[proto.name] = proto;
		});
		//console.log(">>>>extend", con.prototype);
	},
	
/**
Serialize an Array to the callback cb(rec,info) or cb(null,stack) at end given 
a sync/async fetcher( rec, res ).
@param {function} fetch Callback to fetch the data sent to the cb
@param {function} cb Callback to process the fetched data.	
*/
	serial: function (fetch, cb) {
		//Stream( this, {}, (rec, key, res) => {
		this.stream( (rec, key, res) => {
			if ( res )
				fetch( rec, info => {
					cb(rec, info);	// forward results
					res();	// signal record processed w/o stacking any results
				});	

			else 
				cb( null, res );
		});
	},

/**
@param {Function} cb Callback(arg,idx)
@returns this
*/
	any: function ( cb ) {
		var test = false;
		if ( cb )
			if ( typeof cb == "function" )
				this.forEach( val => {
					if ( cb(val) ) test = true;
				});
		
			else
				this.forEach( val => {
					if ( val == cb ) test = true;
				});
				
		else
			this.forEach( val => {
				if ( val ) test = true;
			});
			
		return test;
	},
	
/**
@param {Function} cb Callback(arg,idx)
@returns this
*/
	all: function ( cb ) {
		var test = true;
		if ( cb )
			if ( typeof cb == "function" )
				this.forEach( val => {
					if ( !cb(val) ) test = false;
				});
		
			else
				this.forEach( val => {
					if ( val != cb ) test = false;
				});
				
		else
			this.forEach( val => {
				if ( !val ) test = false;
			});
			
		return test;
	},
	
/**
Index an array using a indexor:

	string of the form "to=from & to=eval & to & ... & !where=eval"
	hash of the form {to: from, ...}
	callback of the form (idx,array) => { ... }

The "!where" clause returns only records having a nonzero eval. 

@param {String|Object|Function} index Indexer 
@param {Object} ctx Context of functions etc 
@returns {Object} Indexed data

@example
[{x:1,y:2},{x:10,y:20}].get("u=x+1&v=sin(y)&!where=x>5",Math)
{ u: [ 11 ], v: [ 0.9129452507276277 ] }

@example
[{x:1,y:2},{x:10,y:20}].get("x")
{ x: [ 1, 10 ] }

@example
[{x:1,y:2},{x:10,y:20}].get("x&mydata=y")
{ mydata: [ 2, 20 ], x: [ 1, 10 ] }

@example
[{x:1,y:2},{x:10,y:20}].get("mydata=[x,y]")
{ mydata: [ [ 1, 2 ], [ 10, 20 ] ] }

@example
[{x:1,y:2},{x:10,y:20}].get("mydata=x+1")
{ mydata: [ 2, 11 ] }

@example
[{x:1,y:2},{x:10,y:20}].get("",{"!all":1})
{ x: [ 1, 10 ], y: [ 2, 20 ] }

@example
[{x:1,y:2},{x:10,y:20}].get("")
[ { x: 1, y: 2 }, { x: 10, y: 20 } ]

@example
[{x:1,y:2},{x:10,y:20}].get("u")
{ u: [ undefined, undefined ] }

@example
[[1,2,3],[10,20,30]].get("1&0")
{ '0': [ 1, 10 ], '1': [ 2, 20 ] }	
*/
	get: function (index,ctx) {
		return getList(this, index, VM.createContext( ctx || {}) )
	},
	
/*
@param {Function} cb Callback(arg,idx)
@returns this
*/
	/*
	put: function ( cb ) {
		const tar = this;
		this.forEach( (arg,i) => tar[i] = cb(arg,i) );
		return this;
	}, */

/*
@param {Function} cb Callback(arg,idx)
@returns this
*/
	/*
	select: function (cb) {
		const rtn = [];
		this.forEach( (arg,i) => {
			if ( res = cb(arg,i) ) rtn.push(res);
		});
		return rtn;
	} */

}, Array.prototype); 

/**
@module ENUMS.String
*/
Copy({
	stream: function (opts,cb) {
		const
			msgs = []; 
		
		this.streamFile( opts, recs => {
			if ( recs )
				recs.forEach( (rec,idx) => {
					cb( rec, idx, msg => {
						if ( msg != undefined ) 
							msgs.push( msg );
					});
				});

			else 	// signal end
				cb(null, msgs);
		});		
	},
	
/**
*/
	replaceSync: function ( pat, cb ) {
		var todo = 0, done=0, rtn = this+"";

		rtn = this.replace( pat, (...args) => {
			const hold = `@@hold${todo++}`;

			//console.log("hold", hold, todo, done);
			cb( args, res => {
				//console.log("res return(", rtn, ")=>", res, done, todo);
				rtn = rtn.replace( hold, res );

				if ( ++done == todo ) 
					cb( rtn ); 
			});
			return hold;
		});

		//console.log("rtn scanned", rtn);
		if ( !todo ) cb( rtn );
	},

/**
Tag url with specified attributes.

@param {String} el tag html element or one of "?&/:=" 
@param {String} at tag attributes = {key: val, ...}
@return {String} tagged results
*/ 
	tag: function (el,at) {
		switch (el) {
			case "/":
			case "?":
			case "&":   // tag a url
				const
					[path,search] = this.split("?"),
					parms = search ? search.split("&") : [],
					keys = {};
				
				var rtn = path;

				parms.forEach( parm => {
					const
						[lhs,rhs] = parm.split("=");
					keys[lhs] = rhs;
				});
				
				if (at)
					Each(Copy(at,keys), (key,val) => {
						rtn += el + key + "=" + val;
						el = "&";
					});
				
				/*
				var rtn = this;

				if (rtn.indexOf("?") >= 0) el = "&";
				
				Each(at, (key,val) => {
					if ( val ) {
						rtn += el + key + "=" + val;
						el = "&";
					}
				});*/

				return rtn;	

			case "[]":
			case "()":
				var rtn = this+el.substr(0,1), sep="";
				
				if (at)
					Each(at, (key,val) => {
						rtn += sep + key + ":" + JSON.stringify(val);
						sep = ",";
					});
				return rtn+el.substr(-1);

			case ":":
			case "=":
				var rtn = this, sep="";
				
				if (at)
					Each(at, (key,val) => {
						rtn += sep + key + el + JSON.stringify(val);
						sep = ",";
					});
				return rtn;

			case "":
				return `<a href="${el}">${this}</a>`;

			default: // tag html

				var rtn = "<"+el+" ";

				if ( at )
					Each( at, (key,val) => {
						rtn += key + "='" + val + "' ";
					});

				switch (el) {
					case "embed":
					case "img":
					case "link":
					case "input":
						return rtn+">" + this;
					default:
						return rtn+">" + this + "</"+el+">";
				}
		}
	},

/**
Parse "$.KEY" || "$[INDEX]" expressions given $ hash.

@param {Object} $ source hash
*/
	parseEval: function ($) {
		try {
			return eval(this+"");
		}

		catch (err) {
Trace("err parse eval", this);
			return err+"";
		}
	},

/**
Run JS against string in specified context.

@param {Object} ctx context hash
*/
	parseJS: function (ctx, cb) {
		try {
			return VM.runInContext( this+"", VM.createContext(ctx || {}));
		}
		catch (err) {
Trace("parseJS", this+"", err, ctx);
			if ( cb ) 
				return cb(err);

			else
				return null;
		}
	},

/**
Return an EMAC "...${...}..." string using supplied context.

@param {Object} query context hash
*/
	parse$: function (ctx) {
		try {
			return VM.runInContext( "`" + this + "`" , VM.createContext(ctx));
		}
		catch (err) {
Trace("parse$ err", err,this);
			return err+"";
		}
	},

/**
Parse string into json or set to default value/callback if invalid json.

@param {Function | Object} def default object or callback that returns default
*/
	parseJSON: function (def) {
		try { 
			return JSON.parse(this);
		}
		catch (err) {  
Trace("parsejson err", this, err);
			return def ? (isFunction(def) ? def(this+"") : def) || null : null;
		}
	},

/**
Parse a "PATH?PARM&PARM&..." url into the specified query, index, flags, or keys hash
as directed by the PARM = ASKEY := REL || REL || _FLAG = VALUE where 
REL = X OP X || X, X = KEY || KEY$[IDX] || KEY$.KEY and returns [path,file,type].

@param {Object} query hash of query keys
@param {Object} index hash of sql-ized indexing keys
@param {Object} flags hash of flag keys
@param {Object} where hash of sql-ized conditional keys
@return {Array} [path,table,type,area,url]
*/
	parsePath: function (query,index,flags,where) { 
		const
			regs = {
				flags: /^\_(.*?)(=)(.*)/,
				index: /(.*?)(:=)(.*)/,
				where: /(.*?)(<=|>=|\!=|\_in=|\_out=)(.*)/,
				other: /(.*?)(<|>|=)(.*)/,

				area: /\/(.*?)\/(.*)/,
				file: /(.*)\.(.*)/
			},
			url = new URL(this, "file:"),
			{pathname,search} = url,
			[x1, area, rem] = pathname.match( regs.area ) || ["", "", pathname.substr(1)],
			[x2, table, type] = rem ? rem.match( regs.file ) || ["", rem, rem.endsWith("/") ? "/" : "" ] : ["","","/"];

		//console.log("parsepath", {search:search, area:area, table: table, type: type, rem: rem});

		["=", "<", "<=", ">", ">=", "!=", "_bin=", "_exp=", "_nlp=","_in=","_out="]
		.forEach( key => where[key] = {} );

		unescape(search).substr(1).split("&").forEach( parm => {
			parm = parm
			.replace( regs.flags, (str,lhs,op,rhs) => {
				try {
					flags[lhs] = JSON.parse(rhs);
				}
				catch (err) {
					flags[lhs] = rhs;
				}
				return "";
			})
			.replace( regs.index, (str,lhs,op,rhs) => {
				index[lhs] = rhs;
				return "";
			})
			.replace( regs.where, (str,lhs,op,rhs) => {
				(where[op] || {})[lhs] = rhs;
				return "";
			})
			.replace( regs.other, (str,lhs,op,rhs) => {
				switch (op) {
					case "=":
						try {
							query[lhs] = JSON.parse(rhs);
						}
						catch (err) {
							query[lhs] = rhs;
						}
						
					default:
						where[op][lhs] = rhs;
				}
				return "";
			});
			
			if ( parm ) 
				index[parm] = "";
		});

		//console.log("parse", [query,index,flags,where]);

		return [pathname,table,type,area,url];
	},

/**
Chunk stream at path by splitting into newline-terminated records.
Callback cb(record) until the limit is reached (until eof when !limit)
with cb(null) at end.

@param {String} path source file
@param {Object} opts {newline,limit} options
@param {Function) cb Callback(record)
*/
	chunkFile: function ( {newline,limit}, cb) {
		const
			path = this+"";

		FS.open( path, "r", (err, fd) => {
			//Trace("chunk", path, err, limit);

			if (err) 
				cb(null);	// signal pipe end

			else {	// start pipe stream
				var 
					pos = 0,
					rem = "",
					run = true,
					src = FS.createReadStream( "", { fd:fd, encoding: "utf8" }),
					sink = new STREAM.Writable({
						objectMode: true,
						write: (bufs,en,sinkcb) => {
							//Trace( "chunk stream", limit,pos,!limit || pos<limit);
							if ( run ) {
								var 
									recs = (rem+bufs).split(newline);

								rem = recs.pop();

								if ( limit )	// chunk limited records
									recs.forEach( (rec,n) => {
										if ( rec ) 
											if ( pos <= limit ) 
												cb(rec);
											
											else {
												//Trace("chunk halting");
												run = false;
											}
										pos++;
									});

								else	// chunk all records
									recs.forEach( (rec,n) => {
										if ( rec ) 
											cb(rec);
										pos++;
									});									
							}
							sinkcb(null);  // signal no errors
						}
					});

				sink
					.on("finish", () => {	// signal complete
						//Trace("chunk finish", limit, pos);
						if (!limit || pos<=limit)	// flush buffer
							if ( rem ) cb(rem);

						cb(null);	// signal complete
					})
					.on("error", err => cb(null) );

				src.pipe(sink);  // start the pipe
			}
		});
	},

/**
Parse a csv/txt/json stream at the specified path dependings on if the
keys is

	[] then record keys are determined by the first header record; 
	[ 'key', 'key', ... ] then header keys were preset; 
	null then raw text records are returned; 
	function then use to parse records.  

The file is chunked using the (newline,limit) chinkFile parameters.  
Callsback cb(record) for each record with cb(null) at end.

@param {String} path source file
@param {Object} opts {keys,comma,newline,limit} options
@param {Function} cb Callback(record || null)
*/
	parseFile: function ( {keys, comma, newline, limit}, cb) {
		const
			path = this+"",
			opts = {newline:newline,limit:limit};

		var 
			pos = 0;

//Trace("parse keys",keys);
		
		if ( keys ) {		// split csv/json file
			const
				parse = keys.forEach 		// define buffer parser
					? (buf,keys) => {		// parse csv buffer
						if ( pos ) {	/// at data row
							var 
								rec = {},
								cols = buf.split(comma);
//Trace("parse cols", cols,keys);
							keys.forEach( (key,m) => rec[key] = cols[m] );
							return rec;
						}

						else {	// at header row so define keys
							if ( buf.charCodeAt(0) > 255 ) buf=buf.substr(1);	// weird
							buf.split(",").forEach( key => keys.push(key) );
//Trace("parse header keys", keys);
							return null;
						}
					}

					: buf => keys( buf );	// parse json buffer

			path.chunkFile( opts, buf => {			// get a record buffer
				if ( buf ) 
					if (rec = parse(buf,keys) ) {	// feed valid record
						pos++;
						cb(rec);
					}
					
					else							// drop invalid record
						pos++;

				else	// forward eof signal 
					cb(null);
			});
		}

		else				// split txt file
			path.chunkFile( opts, buf => {
				if ( buf ) {
					pos++;
					cb(buf);
				}

				else	// forward end signal 
					cb(null);
			});			
	},

	/*
	function streamFile({batch, keys, comma, newline, limit, rekey}, cb) {
		const 
			path = this+"",
			recs = [],
			opts = {keys:keys, limit:limit, comma:comma||",", newline:newline||"\n"},
			rekeyStash = {
				$ctx: {
					$tag: (str,el,tags) => str.tag(el,tags),
					$link: (str,ref) => str.tag("a",{href:ref}),
					$grab: str => str.replace(/.*\#(.*)\#.* /, (str,arg) => arg )  <<<* /
					//+ $.get(list,"[from || regexp || (JS) || !rem] => [to || !test] || INDEX, ..."  )  
					//+ $.resize(256,256)
					//+ $.greyscale()
					//+ $.sym( {maps:{x:'RGB',y:'L'}, limits:{cols:N,rows:N}, levels:{x:N,y:N}} )					
				}
			};

		var
			pos = 0;

		path.parseFile( opts, rec => {
			if ( rec ) 
				if ( !batch || recs.length<batch ) 
					recs.push( rec );

				else {	// flush batch
					pos += recs.length;

					const 
						recBatch = recs.get(rekey,rekeyStash);

					if ( recBatch.length ) cb( recBatch );
					recs.length = 0;
				}

			else { // flush batch then forward end signal
				//Trace("!!!!!!!!!!!!eos", recs.length);
				pos += recs.length;

				const 
					recBatch = recs.get(rekey,rekeyStash);

				if ( recBatch.length ) cb( recBatch );
				//cb(recs.get(rekey,rekeyStash), pos+=recs.length); // flush batch
				//Trace("!!!!!!!!!!!!eos end");
				cb(null);	// signal end
			}
		});
	},*/

/**
Stream file at path containing comma delimited values.  The file is split using the (keys,comma) 
file splitting parameters, and chunked using the (newline,comma) file chunking parameters. Callsback 
cb( [record,...] ) with the record batch or cb( null ) at end.  

@param {String} path source file
@param {Object} opts {keys,comma,newline,limit,batch} options
@param {Function} cb Callback( [record,...] || null )
*/
	streamFile: function ( {batch, keys, comma, newline, limit}, cb) {
		const 
			path = this+"",
			recs = [],
			opts = {keys:keys, limit:limit, comma:comma||",", newline:newline||"\n"};		
		var
			pos = 0;

		path.parseFile( opts, rec => {
			if ( rec ) 
				if ( !batch || recs.length<batch ) 
					recs.push( rec );

				else {	// flush batch
					pos += recs.length;
					if ( recs.length ) cb( recs );
					recs.length = 0;
				}

			else { // flush batch then forward end signal
				pos += recs.length;
				if ( recs.length ) cb( recs );
				cb(null);	// signal end
			}
		});
	},

/**
Trace message to console with optional request to place into syslogs
@param {String} msg message to trace
@param {Object} req request { sql, query, client, action, table } 
@param {Function} res response callback(msg)
*/
	trace: function (...args) {	
		function cpu() {
			var sum = 0, util = OS.loadavg();
			for ( var n=0, N = util.length; n<N; n++) sum += util[0];
			return sum / N;
		}

		function mem() {
			return OS.freemem() / OS.totalmem();
		}

		const
			node = CLUSTER.isMaster ? 0 : CLUSTER.worker.id;

		if ( this.endsWith("*") ) {
			const 
				[req] = args || [],
				{ sql, query, client, action, table } = req || {};

			if ( sql ) 
				sql.query( "INSERT INTO openv.syslogs SET ? ", {
					Node: OS.hostname()+"."+node,
					Client: client,
					Table: table,
					At: new Date(),
					Case: query.Name || "",
					Action: action,
					Module: this,
					cpuUtil: cpu(),
					memUtil: mem(),
					Message: this
				});
		}
		
		else
			console.log( node + ">>>"  + this.toUpperCase(), args[0] );		
	},

/**
Serialize this String to the callback(results) given a sync/asyn fetcher(rec,res) where
rec = {ID, arg0, arg1, ...} contains args produced by regex.  Provide a unique placeholder
key to back-substitute results.

@example

	"junkabc;junkdef;"
	.serial( (rec,cb) => cb("$"), /junk([^;]*);/g, "@tag", msg => console.log(msg) )

produces:

	"$$"
*/
	serial: function ( fetcher, regex, key, cb ) {  //< callback cb(str) after replacing regex using fetcher( rec, (ex) => "replace" ) and string place holder key
		var 
			recs = [],
			results = this.replace( regex, (arg0, arg1, arg2, arg3, arg4) => {  // put in place-holders
				recs.push( {ID: recs.length, arg0:arg0, arg1:arg1, arg2:arg2, arg3:arg3, arg4:arg4} );
				return key+(recs.length-1);
			});

		recs.serial( fetcher, (rec,info) => {  // update place-holders with info 
			if (rec) 
				results = results.replace( key+rec.ID, str => info );

			else 
				cb( results );
		});

	},

/**
Fetch using supplied url.
*/
	get: function (cb) {
		Fetch(this,cb);
	}
		
}, String.prototype);

/**
Create a clock object with specified trace switch, every interval, on-off times, and 
start date.  See the clock tick method for more information.

@module ENUMS.Clock
@param {Boolean} trace tracking switch
@param {String | Float} every tick interval
@param {Float} on on-time or 0 implies infinity
@param {Float} off off-time or 0
@param {Date} start start date
*/
function Clock(trace,every,on,off,start) {
	this.next = start ? new Date( start ) : new Date();
	this.on = on || Infinity;
	this.off = off || 0;
	this.epoch = this.off;
	this.every = every.replace( /[0-9.]/g, "");		// N [sec|min|hour|....]
	this.trace = trace;
	this.cycle = off ? on+off-2 : 1; //(on>=2) ? on+off-2 : 0;  // Internal clock cycle = CYCLE - 2 = ON+OFF-2.
	this.step = parseFloat(every) || 1;
	this.util = on ? on/(on+off) : 0;
	
	if ( trace ) Trace(this);
}

Copy({
/**
@returns {Date} Current clock time
*/	
	now: function () {
		return new Date( this.next );
	},
	
/**
Return the wait time to next event, with callback(wait,next) when at snapshot events.

Example below for ON = 4 and OFF = 3 steps of length clock.every = sec|min|hour|...

Here S|B|* indicates the end of snapshot|batch|start events.  The clock starts on epoch = OFF 
with a wait = 0.  The clock's host has 1 step to complete its batch tasks, and OFF steps to 
complete its snapshot tasks.  Here, the work CYCLE = ON+OFF with a utilization = ON/CYCLE.  
Use OFF = 0 to create a continious process.  

		S			*	B	B	B	S			*	B	B	B
		|			|	|	|	|	|			|	|	|	|		...
		|			|	|	|	|	|			|	|	|	|		
		x-->x-->x-->x-->x-->x-->x-->x-->x-->x-->x-->x-->x-->x-->x-->x-->x-->x
epoch	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17

		|<-- OFF -->|<---- ON ----->|

@param {Function} cb Callback
@returns {Date} Wait time
*/	
	tick: function (cb) {		// advance clock, return waitTime, and callback cb(nextTime) on state-change epochs 
		
		const
			{ max,trunc } = Math,
			{ step,every,on,off,epoch,trace,next,start,cycle } = this,
			now = new Date(),
			wait = max( 1e3, next.getTime() - now.getTime() );

		if (trace)
			Trace( `e=${epoch} s=${step} w=${trunc(wait)} n=${next}` );

		if ( cycle > 1 )	// using on-off process
			if ( epoch % cycle == 0 ) {
				this.step = off;
				if (trace)
					Trace( `e=${epoch} s=${step} w=${trunc(wait)} n=${next}` );
				if (cb) cb( wait, next );
			}
		
			else			// continious process
				this.step = 1;
		
		this.epoch += this.step;
		
		switch ( every ) {		// advance next job epoch
			case "yr":
			case "year": 	next.setFullYear( next.getFullYear() + this.step ); break;
			case "wk":
			case "week": 	next.setDate( next.getDate() + 7*this.step ); break;
			case "day": 	next.setDate( next.getDate() + this.step ); break;
			case "hr":
			case "hour": 	next.setHours( next.getHours() + this.step ); break;
			case "min":
			case "minute": 	next.setMinutes( next.getMinutes() + this.step ); break;
			case "sec":
			case "second": 	next.setSeconds( next.getSeconds() + this.step ); break;

			case "monday": 
			case "tuesday":
			case "wednesday":
			case "thursday":
			case "friday":
			case "saturday":
			case "sunday":
			case "mon":
			case "tue":
			case "wed":
			case "thr":
			case "fri":
			case "sat":
			case "sun":
				next.setDate( next.getDate() + 7 - next.getDay() + 1); 
				this.every = "week";	// now switch to weeks
				// Trace("clock next", this.next, this.every);
				break;
				
			default:
				next.setTime( next.getTime() + every*1e3 );

		}

		return wait;
	}
}, Clock.prototype);

// Establish connections to services

var
	mysqlThreads = 0;

const
	{mysqlCon, neo4jCon, txmailCon, rxmailCon} = Copy({
		mysqlCon: mysqlLogin 
					? MYSQL.createPool( Copy({
						host: mysqlLogin.hostname,
						user: mysqlLogin.username,
						password: mysqlLogin.password,
						//port: mysqlLogin.port
					}, mysqlOpts) )
					// reserved to test thread depth to protect against DOS attacks
					.on("acquire", () => mysqlThreads ++ )
					.on("release", () => mysqlThreads -- )
		
					: null,
		
		neo4jCon: neo4jLogin 
					? NEO4J.driver( "bolt://localhost", NEO4J.auth.basic( 
						neo4jLogin.username, neo4jLogin.password ), neo4jOpts )
		
					: null,
		
		rxmailCon: rxmailLogin
					? new IMAP( Copy({
						user: rxmailLogin.username,
						password: rxmailLogin.password,
						host: rxmailLogin.hostname,
						port: rxmailLogin.port
					}, rxmailOpts) )
		
					: null,

		txmailCon: txmailLogin 
					? MAIL.createTransport({
						host: txmailLogin.hostname,
						port: txmailLogin.port,
						auth: {
							user: txmailLogin.username,
							pass: txmailLogin.password
						}
					})

					: {
						sendMail: (opts, cb) => {
							Trace("email", opts);   // -r "${opts.from}" 
			
							if ( false ) //DEBE.watchMail
								sqlThread( sql => {
									sql.query("INSERT INTO openv.email SET ?", {
										To: opts.to,
										Body: opts.body,
										Subject: opts.subject,
										Send: false,
										Remove: false
									});
								});

							else
								exec(`echo -e "${opts.body}\n" | mailx -s "${opts.subject}" ${opts.to}`, err => {
									cb( err );
									//Trace("MAIL "+ (err || opts.to) );
								});
						}
					}

	}, ENUMS);

if (rxmailCon)					// Establish server's email inbox	
	rxmailCon.connect( err => {  // login cb
		if (err) Trace(err);

		rxmailCon.openBox('INBOX', true, (err,mailbox) => {

			if (err) Trace(err);

			rxmailCon.search([ 'UNSEEN', ['SINCE', 'May 20, 2012'] ], (err, results) => {

				if (err) Trace(err);

				rxmailCon.Fetch(results, { 
					headers: ['from', 'to', 'subject', 'date'],
					cb: fetch => {
						fetch.on('message', msg => {
							Trace('Saw message no. ' + msg.seqno);
							msg.on('headers', hdrs => {
								Trace('Headers for no. ' + msg.seqno + ': ' + hdrs);
							});
							msg.on('end', () => {
								Trace('Finished message no. ' + msg.seqno);
							});
						});
					}
				}, err => {
					if (err) throw err;
					Trace('Done fetching all messages!');
					rxmailCon.logout();
				});
			});
		});
	});

function NEOCON(trace) {
	this.trace = trace || false;
}

const
	{neoThread, sqlThread} = Copy({
		neoThread: cb => {
			if ( cb )
				return neo4jCon ? cb( new NEOCON( false ) ) : null;

			else
				return neo4jCon ? true : false;
		},

		// callback cb(sql) with a sql connection
		sqlThread: cb => {	
			
			//console.log("===========>thread", mysqlCon.constructor.prototype);
			
			if ( cb )
				return mysqlCon ? cb( mysqlCon ) : null;

			else
				return mysqlCon ? true : false;
		}
	}, ENUMS);

//================== Unit testing

async function LexisNexisTest(N,endpt,R,cb) {
	const start = new Date(), {random,trunc} = Math;
	var done = 0;
	Trace(start);
	for ( var n=0; n<N; n++) 
		Fetch(endpt + (R?trunc(random()*R):""), res => {
			//Trace(n,done,N);
			if ( ++done == N ) {
				var 
					stop = new Date(),
					mins = (stop-start)/1e3/60,
					rate = N/mins;
				
				Trace(stop, mins, "mins", rate, "searches/min");
				
				if (cb) cb(res);
			}
		});
}

Start("enums", Copy({
	"??": () => 
		Trace("", JSON.stringify({
			logins: {
				mysql: mysqlLogin,
				neo4j: neo4jLogin,
				txmail: txmailLogin,
				rxmail: rxmailLogin
			},
			connections: {
				mysql: mysqlCon ? true : false,
				neo4j: neo4jCon  ? true : false,
				txmail: txmailCon ? true : false,
				rxmail: rxmailCon ? true : false
			}
		})),
	
	EMT: () => 
	 	ENUMS.txmailCon.sendMail({from: "noreply@totem.org", to: "brian.d.james@comcast.net", subject: "this is a test", body: "hello"}),
	
	ECOPY: () =>
	  	Trace("", {
			shallowCopy: Copy( {a:1,b:2}, {} ),
			deepCopy: Copy({ a:1,"b.x":2 }, {b:{x:100}}, ".")
		}),
			
	ELN1: () =>
		LexisNexisTest(1e3, 'lex://services-api.lexisnexis.com/v1/News?$search=rudolph'),
		
	ELN2: () =>
		LexisNexisTest(1e3, 'lex://services-api.lexisnexis.com/v1/News?$search=rudolph&$expand=Document'),
		
	ELN3: () =>
		LexisNexisTest(1e3, 'lex://services-api.lexisnexis.com/v1/News?$search=rudolph'),
		
	ELN4: () =>
		LexisNexisTest(1, 'lex://services-api.lexisnexis.com/v1/News?$search=rudolph', 0, res => {
			//Trace("res=>", res);
			var r = JSON.parse(res);
			//Trace( Object.keys(r) );
			if ( rec = r.value[0] ) 
				Trace( "fetch docendpt", rec['DocumentContent@odata.mediaReadLink'] , Object.keys(rec), rec.Overview, rec.Date );
			
			Fetch( 'lex://services-api.lexisnexis.com/v1/' + rec['DocumentContent@odata.mediaReadLink'] , doc => {
				//    'lexis://services-api.lexisnexis.com/v1/MEDIALINK
				Trace( "doc", doc );
			});
			
		}),
		
	ELN5: () =>
		Fetch( 'lexis://services-api.lexisnexis.com/v1/News?$search=rudolph', doc => {
			Trace( "doc", doc );
		}),
		
	EX1: () => 
		"testabc; testdef;".replaceSync(/test(.*?);/g, (args,cb) => {
			
			if ( cb ) {
				console.log("args", args);
				cb( "#"+args[1] );
			}
			
			else
				console.log("final", args);
		}),
		
	EX2: () =>
		"testabc; testdef;".replaceSync(/test(.*?);/g, (args,cb) => {
			//console.log(args);
			
			if ( cb ) 
				Fetch( "http://localhost:8080/nets.txt?x:=Name", txt => {
					//console.log("fetch", args, txt);
					cb( "#"+args[1] );
					//cb( "#"+txt );
				});
			
			else
				console.log("final", args); 
		})
}, (cmd,ctx) => {
	VM.runInContext(cmd, VM.createContext(ctx));
}));

// UNCLASSIFIED
