/**********************************
 *** Common routines from k8.js ***
 **********************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/*************/

function fm_log2tbl(args)
{
	if (args.length == 0) {
		warn("Usage: k8 fermi2.js log2tbl <file1.log> [...]");
		exit(1);
	}
	var buf = new Bytes();
	for (var i = 0; i < args.length; ++i) {
		var fn = args[i];
		var f = new File(fn);
		var prev = 0, curr = 0;
		while (f.readline(buf) >= 0) {
			var m, line = buf.toString();
			if ((m = /mr_restore.*\((\d+)/.exec(line)) != null) {
				prev = m[1];
			} else if ((m = /main_ropebwt2.*\((\d+)/.exec(line)) != null) {
				curr = m[1];
			}
		}
		f.close();
		print(fn, curr, prev);
	}
}

function fm_id2sam(args)
{
	if (args.length == 0) {
		warn("Usage: k8 fermi2.js id2sam <sample.tbl> <matches.txt>");
		exit(1);
	}
	var buf = new Bytes();
	var f = new File(args[0]);
	var sample = [];
	while (f.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		sample.push([parseInt(t[1]), t[0]]);
	}
	f.close();
	sample.sort(function(a,b) {return a[0]-b[0];})
	f = args.length > 1? new File(args[1]) : new File();
	while (f.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0] == 'EM') {
			var ids = [];
			for (var i = 4; i < t.length; ++i) {
				var s = t[i].split(":");
				var id = parseInt(s[0]);
				if (id >= sample[sample.length-1][0])
					throw Error("Invalid sequence ID "+id);
				ids.push(id);
			}
			ids.sort(function(a,b){return a-b;});
			for (var i = 0; i < ids.length; ++i) {
				var id = ids[i];
				var start = 0, end = sample.length;
				while (start < end) {
					var mid = (start + end) >> 1;
					if (sample[mid][0] == id) start = end = mid - 1;
					else if (sample[mid][0] < id) start = mid + 1;
					else end = mid;
				}
				t[i+4] = sample[start][1];
			}
			print(t.join("\t"));
		} else print(buf);
	}
	f.close();
}

function fm_cnt2hist(args)
{
	var buf = new Bytes();
	var f = args.length? new File(args[0]) : new File();
	var tbl = [];
	while (f.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var x = parseInt(t[1]);
		if (tbl[x] == null) tbl[x] = 0;
		++tbl[x];
	}
	f.close();
	buf.destroy();
	for (var i = 0; i < tbl.length; ++i)
		if (tbl[i] != null) print(i, tbl[i]);
}

function fm_unreg(args)
{
	var c, min_len = 61, min_utg_len = 0, gap = 1;
	while ((c = getopt(args, 'g:l:u:')) != null) {
		if (c == 'g') gap = parseInt(getopt.arg);
		else if (c == 'l') min_len = parseInt(getopt.arg);
		else if (c == 'u') min_utg_len = parseInt(getopt.arg);
	}

	var f = args.length == getopt.ind? new File() : new File(args[getopt.ind]);
	var b = new Bytes();

	function process(name, len, a)
	{
		var s = 0, e = -1, c = [];
		if (len < min_utg_len) return;
		for (var i = 0; i < a.length; ++i) {
			if (a[i][0] > e) {
				if (e > 0) c.push([s, e]);
				s = a[i][0]; e = a[i][1];
			} else if (a[i][1] > e) e = a[i][1];
		}
		if (e > 0) c.push([s, e]);
		s = 0;
		for (var i = 0; i < c.length; ++i) {
			if (c[i][0] > s && s >= min_len) print(name, s, c[i][0]);
			s = c[i][1];
		}
		if (len > s && len - s >= min_len) print(name, s, len);
		c = [];
	}

	var name, len, a;
	while (f.readline(b) >= 0) {
		var m, line = b.toString();
		if ((m = /^SQ\s(\S+)\s(\d+)/.exec(line)) != null) {
			name = m[1]; len = parseInt(m[2]); a = [];
		} else if ((m = /^EM\s(\d+)\s(\d+)/.exec(line)) != null) {
			var st = parseInt(m[1]), en = parseInt(m[2]);
			if (en - st >= min_len) a.push([st + gap, en - gap]);
		} else if (/^\/\//.test(line)) {
			process(name, len, a);
	}
	}

	b.destroy();
	f.close();
}

/***********************
 *** Main() function ***
 ***********************/

function main(args)
{
	if (args.length == 0) {
		print("\nUsage:    k8 fermi2.js <command> [arguments]\n");
		print("Commands: log2tbl   extract sample info from ropebwt2 log");
		print("          id2sam    relate sequence index to sample info");
		print("          cnt2hist  k-mer count histogram");
		print("          unreg     identify unaligned regions from match output");
		print("");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'log2tbl') fm_log2tbl(args);
	else if (cmd == 'id2sam') fm_id2sam(args);
	else if (cmd == 'cnt2hist') fm_cnt2hist(args);
	else if (cmd == 'unreg') fm_unreg(args);
	else warn("Unrecognized command");
}

main(arguments);
