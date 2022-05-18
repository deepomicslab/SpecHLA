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

function intv_ovlp3(intv, bits, to_merge)
{
	if (typeof bits == "undefined") bits = 13;
	if (typeof to_merge == "undefined") to_merge = true;
	intv.sort(function(a,b) {return a[0]-b[0];});
	if (to_merge) { // merge overlapping regions
		var j = 0;
		for (var i = 1; i < intv.length; ++i) {
			if (intv[j][1] > intv[i][0])
				intv[j][1] = intv[j][1] > intv[i][1]? intv[j][1] : intv[i][1];
			else intv[++j] = intv[i].slice(0);
		}
		intv.length = j + 1;
	}
	//
	var sum = 0;
	for (var i = 0; i < intv.length; ++i)
		sum += intv[i][1] - intv[i][0];
	// create the index
	var idx = [], max = 0;
	for (var i = 0; i < intv.length; ++i) {
		var b = intv[i][0]>>bits;
		var e = (intv[i][1]-1)>>bits;
		if (b != e) {
			for (var j = b; j <= e; ++j)
				if (idx[j] == null) idx[j] = i;
		} else if (idx[b] == null) idx[b] = i;
		max = max > e? max : e;
	}
	return [function(_b, _e) { // closure
		var x = _b >> bits;
		if (x > max) return null;
		var off = idx[x];
		if (off == null) {
			var i;
			for (i = ((_e - 1) >> bits) - 1; i >= 0; --i)
				if (idx[i] != null) break;
			off = i < 0? 0 : idx[i];
		}
		for (var i = off; i < intv.length && intv[i][0] < _e; ++i)
			if (intv[i][1] > _b) return intv[i];
		return null;
	}, sum]
}

function intv_ovlp2(intv, bits)
{
	return intv_ovlp3(intv, bits, true);
}

function intv_ovlp(intv, bits)
{
	return intv_ovlp2(intv, bits)[0];
}

Math.lgamma = function(z) {
	var x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return Math.log(x) - 5.58106146679532777 - z + (z-0.5) * Math.log(z+6.5);
}

Math.fisher_exact = function(n11, n12, n21, n22)
{
	function lbinom(n, k) {
		if (k == 0 || n == k) return 0;
		return Math.lgamma(n+1) - Math.lgamma(k+1) - Math.lgamma(n-k+1);
	}

	function hypergeo(n11, n1_, n_1, n) {
		return Math.exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
	}

	function hypergeo_acc(n11, n1_, n_1, n, aux) {
		if (n1_ || n_1 || n) {
			aux.n11 = n11; aux.n1_ = n1_; aux.n_1 = n_1; aux.n = n;
		} else { // then only n11 changed; the rest fixed
			if (n11%11 && n11 + aux.n - aux.n1_ - aux.n_1) {
				if (n11 == aux.n11 + 1) { // incremental
					aux.p *= (aux.n1_ - aux.n11) / n11
						* (aux.n_1 - aux.n11) / (n11 + aux.n - aux.n1_ - aux.n_1);
					aux.n11 = n11;
					return aux.p;
				}
				if (n11 == aux.n11 - 1) { // incremental
					aux.p *= aux.n11 / (aux.n1_ - n11)
						* (aux.n11 + aux.n - aux.n1_ - aux.n_1) / (aux.n_1 - n11);
					aux.n11 = n11;
					return aux.p;
				}
			}
			aux.n11 = n11;
		}
		aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n);
		return aux.p;
	}

	var i, j, max, min;
	var p, q, left, right, two;
	var _aux = { n11:0, n1_:0, n_1:0, n:0, p:0. };
	var n1_, n_1, n;

	n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
	max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
	min = n1_ + n_1 - n;
	if (min < 0) min = 0; // min n11, for left tail
	if (min == max) return [1., 1., 1.]; // no need to do test
	q = hypergeo_acc(n11, n1_, n_1, n, _aux); // the probability of the current table
	// left tail
	p = hypergeo_acc(min, 0, 0, 0, _aux);
	for (left = 0., i = min + 1; p < 0.99999999 * q; ++i) // loop until underflow
		left += p, p = hypergeo_acc(i, 0, 0, 0, _aux);
	--i;
	if (p < 1.00000001 * q) left += p;
	else --i;
	// right tail
	p = hypergeo_acc(max, 0, 0, 0, _aux);
	for (right = 0., j = max - 1; p < 0.99999999 * q; --j) // loop until underflow
		right += p, p = hypergeo_acc(j, 0, 0, 0, _aux);
	++j;
	if (p < 1.00000001 * q) right += p;
	else ++j;
	// two-tail
	two = left + right;
	if (two > 1.) two = 1.;
	// adjust left and right
	if (Math.abs(i - n11) < Math.abs(j - n11)) right = 1. - left + q;
	else left = 1.0 - right + q;
	return [two, left, right];
}

/*****************************
 *** Compare two BED files ***
 *****************************/

function read_bed(fn, shift)
{
	var f = fn == '-'? new File() : new File(fn);
	var b = new Bytes();
	var reg = {}, idx = {};
	while (f.readline(b) >= 0) {
		var t = b.toString().split("\t");
		if (reg[t[0]] == null) reg[t[0]] = [];
		reg[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
	}
	var n_rec = 0, sum = 0;
	for (var chr in reg) {
		var a = intv_ovlp2(reg[chr], shift)
		idx[chr] = a[0], sum += a[1];
		n_rec += reg[chr].length;
	}
	b.destroy();
	f.close();
	return [idx, n_rec, sum];
}

function b8_bedovlp(args)
{
	if (args.length < 2) {
		print("Usage: k8 hapdip.js bedovlp <read.bed> <inRAM.bed>");
		exit(1);
	}
	var n_rec;
	var f1 = new File(args[0]);
	var tmp = read_bed(args[1]);
	var bed = tmp[0], n_rec = tmp[1];
	var b1 = new Bytes();
	var c = [0, 0, 0];
	while (f1.readline(b1) >= 0) {
		var t = b1.toString().split("\t", 3);
		t[1] = parseInt(t[1]); t[2] = parseInt(t[2]);
		if (bed[t[0]] == null) ++c[0];
		else if (bed[t[0]](t[1], t[2])) ++c[1];
		else ++c[2];
	}
	b1.destroy();
	f1.close();
	print(args[0], args[1], c[0], c[1], c[2], n_rec);
}

/******************************
 * Compare multiple BED files *
 ******************************/

function b8_bedcmpm(args)
{
	function cmpfunc(a, b) {
		if (a[0] < b[0]) return -1;
		if (a[0] > b[0]) return 1;
		return a[1] - b[1];
	}

	var c, win_size = 0, print_joint = false;
	while ((c = getopt(args, 'pw:')) != null)
		if (c == 'w') win_size = parseInt(getopt.arg);
		else if (c == 'p') print_joint = true;

	var bed = [];
	var buf = new Bytes();
	var n = args.length - getopt.ind;
	for (var i = getopt.ind; i < args.length; ++i) {
		var file = new File(args[i]);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			t[1] = parseInt(t[1]);
			if (t.length > 2) t[2] = parseInt(t[2]);
			else t[2] = t[1], --t[1];
			t[1] = t[1] < win_size? 0 : t[1] - win_size;
			t[2] += win_size;
			bed.push([t[0], t[1], t[2], 1 << (i-getopt.ind)]);
		}
		file.close();
	}
	bed.sort(cmpfunc);
	var labels = [];
	for (var i = 0; i < 1<<n; ++i) labels[i] = 0;
	var last = null, start = 0, end = 0, label = 0;
	for (var j = 0; j < bed.length; ++j) {
		if (bed[j][0] != last || bed[j][1] >= end) { // no overlap
			if (last != null) {
				if (print_joint) print(last, start, end, label);
				++labels[label];
			}
			last = bed[j][0], start = bed[j][1], end = bed[j][2], label = bed[j][3];
		} else end = end > bed[j][2]? end : bed[j][2], label |= bed[j][3];
	}
	if (print_joint) print(last, start, end, label);
	++labels[label];
	for (var i = 0; i < 1<<n; ++i)
		print(i, labels[i]);
}

/******************************************
 *** Circular Binary Segmentation (CBS) ***
 ******************************************/

function cnv_cbs1(a0, start, end, winsize) // find the best cut for a single segment
{
	var a = a0.slice(start, end);
	a.sort(function(x,y){return x[1]-y[1]}); // sort by value
	var same = 1, T = 0;
	for (var i = 1; i <= a.length; ++i) { // compute the rank
		if (i == a.length || a[i-1][1] != a[i][1]) {
			var rank = i - .5 * (same - 1);
			for (var j = i - same; j <= i - 1; ++j)
				a[j][2] = rank;
			if (same > 1) T += (same - 1.) * same * (same + 1.), same = 1; // T to correct for ties
		} else ++same;
	}
	var z = ((a.length + 1) - T / a.length / (a.length - 1)) / 12.;
	a.sort(function(x,y){return x[0]-y[0]});
	var s = [], sum = 0;
	s[0] = 0;
	for (var i = 0; i < a.length; ++i) s[i+1] = s[i] + a[i][2];
	var max = [0, -1, -1, 0];
	for (var j = 1; j <= a.length; ++j) {
		for (var i = j > winsize? j - winsize : 0; i < j; ++i) { // compute the standardized Mann-Whitney U
			var n = j - i, m = a.length - n;
			var x = (s[j] - s[i] - .5 * n * (a.length + 1)) / Math.sqrt(m * n * z);
			var y = x > 0? x : -x;
			if (y > max[0]) max = [y, i, j, x]; // this is the best cut
		}
	}
	return max;
}

function cnv_cbs(chr, a, thres, winsize) // recursively cut
{ // FIXME: it does not consider the edge effect; should not be a big issue
	function print_intv(b, e) { // print an interval
		if (b == e) return;
		var t = a.slice(b, e);
		t.sort(function(x,y){return x[1]-y[1]}); // sort to compute the median
		print(chr, a[b][0] - 1, a[e-1][0], t.length, t[Math.floor(t.length/2)][1]);
	}

	warn("Processing chromsome " + chr + "...");
	a.sort(function(x,y){return x[0]-y[0]});
	var queue = [[0, a.length, -1]];
	while (queue.length) {
		var e = queue.shift();
		var start = e[0], end = e[1];
		var cut = cnv_cbs1(a, start, end, winsize);
		if (cut[0] < -thres || cut[0] > thres) { // a significant cut; then cut further
			if (cut[1] >= 4) queue.push([start, start + cut[1]]);
			else print_intv(start, start + cut[1]);
			if (cut[2] - cut[1] >= 4) queue.push([start + cut[1], start + cut[2]]);
			else print_intv(start + cut[1], start + cut[2]);
			if (end - start - cut[2] >= 4) queue.push([start + cut[2], end]);
			else print_intv(start + cut[2], end);
		} else print_intv(start, end); // no significant cut can be found
	}
}

function b8_cbs(args)
{
	var c, winsize = 200, thres = 3.891; // 3.891 is equivalent to P-value 1e-4
	while ((c = getopt(args, 't:w:')) != null) {
		if (c == 't') thres = parseFloat(getopt.arg);
		else if (c == 'w') winsize = parseInt(getopt.arg);
	}

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var last_chr = null, a = [];

	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0] != last_chr) {
			if (a.length) cnv_cbs(last_chr, a, thres, winsize);
			last_chr = t[0];
			a = [];
		}
		a.push([parseInt(t[1]), parseFloat(t[2])]);
	}
	if (a.length) cnv_cbs(last_chr, a, thres, winsize);

	buf.destroy();
	file.close();
}

/****************************
 *** Parse one-sample VCF ***
 ****************************/

function b8_allele_start(ref, alt)
{
	var en_l, st;
	for (en_l = 0; en_l < alt.length && en_l < ref.length; ++en_l)
		if (alt.charAt(alt.length - 1 - en_l) != ref.charAt(ref.length - 1 - en_l))
			break;
	for (st = 0; st < alt.length - en_l && st < ref.length - en_l; ++st) // skip matching bases at the left hand
		if (alt.charAt(st) != ref.charAt(st))
			break;
	return st;
}

function b8_parse_vcf1(t) // t = vcf_line.split("\t")
{
	if (t.length < 6) return null;
	var a = [];
	t[1] = parseInt(t[1]) - 1; t[3] = t[3].toUpperCase(); t[4] = t[4].toUpperCase(); t[5] = parseFloat(t[5]);
	var s = t[4].split(","); // list of ALT alleles
	// find the genotype
	var gt, match = /^(\d+)(\/|\|)(\d+)/.exec(t[9]);
	if (match == null) { // special-casing samtools, as for single-sample, it effectively gives the genotype at FQ
		var m2 = /FQ=([^;\t]+)/.exec(t[7]);
		if (m2 == null) {
			m2 = /^(\.|\d+)(\/|\|)(\.|\d+)/.exec(t[9]);
			if (m2 != null) gt = m2[1] == m2[3]? 1 : 0; // NOTE: ./. is treated as a hom
			else gt = -1;
		} else gt = parseFloat(m2[1]) > 0? 1 : 0;
	} else {
		gt = parseInt(match[1]) != parseInt(match[3])? 1 : 0;
		if (match[1] == 0 && match[2] == 0) return [];
	}
	// get CIGAR for freebayes
	var m3 = /CIGAR=([^;\t]+)/.exec(t[7]);
	var cigar = m3 != null? m3[1].split(",") : [];
	if (cigar.length && cigar.length != s.length) throw Error("Inconsistent ALT and CIGAR");
	// loop through each ALT allele
	for (var i = 0; i < s.length; ++i) {
		var type; // 0=ts, 1=tv, 2=mnp, 3=ins, 4=del
		if (/^<\S+>/.test(s[i])) continue;
		if (t[3].length == 1 && s[i].length == 1) { // SNP; this is a special case of the last "else", but is faster
			type = ((t[3] == 'A' && s[i] == 'G') || (t[3] == 'G' && s[i] == 'A') || (t[3] == 'C' && s[i] == 'T') || (t[3] == 'T' && s[i] == 'C'))? 0 : 1;
			a.push([t[5], type, gt, t[1], 0]);
		} else if (cigar.length) { // MNP or INDEL from freebayes
			var x = 0, y = 0;
			var m4, re = /(\d+)([MXID=])/g;
			while ((m4 = re.exec(cigar[i])) != null) {
				var l = parseInt(m4[1]);
				if (m4[2] == 'X') {
					for (var j = 0; j < l; ++j) {
						u = t[3].charAt(x+j), v = s[i].charAt(y+j);
						type = ((u == 'A' && v == 'G') || (u == 'G' && v == 'A') || (u == 'C' && v == 'T') || (u == 'T' && v == 'C'))? 0 : 1;
						a.push([t[5], type, gt, t[1] + x, 0]);
					}
					x += l, y += l;
				} else if (m4[2] == 'I') {
					a.push([t[5], 3, gt, t[1] + x,  l, s[i].substr(y, l)]);
					y += l;
				} else if (m4[2] == 'D') {
					a.push([t[5], 4, gt, t[1] + x, -l, t[3].substr(x, l)]);
					x += l;
				} else if (m4[2] == 'M' || m4[2] == '=') x += l, y += l;
			}
		} else { // MNP or INDEL from Platypus and others
			var st = b8_allele_start(t[3], s[i]);
			var d = s[i].length - t[3].length, rs, as;
			if (d > 0) { // insertion
				a.push([t[5], 3, gt, t[1] + st, d, s[i].substr(st, d)]);
				rs = st, as = st + d;
			} else if (d < 0) { // deletion
				a.push([t[5], 4, gt, t[1] + st, d, t[3].substr(st, -d)]);
				rs = st + (-d), as = st;
			} else rs = as = st;
			for (var j = 0; j < t[3].length - rs; ++j) { // check the rest
				var u = t[3].charAt(rs + j), v = s[i].charAt(as + j);
				if (u != v) {
					type = ((u == 'A' && v == 'G') || (u == 'G' && v == 'A') || (u == 'C' && v == 'T') || (u == 'T' && v == 'C'))? 0 : 1;
					a.push([t[5], type, gt, t[1] + j + rs, 0]);
				}
			}
		}
	}
	return a; // [qual, ts/tv/ins/del, gt, pos, indelLen, indelSeq]
}

/*****************
 *** VCF stats ***
 *****************/

function b8_qst1(args)
{
	var c, b = 0.02, show_hdr = false, min_q = 0, indel_bed = false, snp_pos = false;
	while ((c = getopt(args, "SGHb:q:")) != null)
		if (c == 'b') b = parseFloat(getopt.arg);
		else if (c == 'q') min_q = parseFloat(getopt.arg);
		else if (c == 'H') show_hdr = true;
		else if (c == 'G') indel_bed = true;
		else if (c == 'S') snp_pos = true;

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var a = [];
	while (file.readline(buf) >= 0) {
		if (buf[0] == 35) continue;
		var t = buf.toString().split("\t");
		var u = b8_parse_vcf1(t);
		if (u == null || u.length == 0 || u[0][0] < min_q) continue;
		if (snp_pos) {
			for (var i = 0; i < u.length; ++i)
				if (u[i][1] == 0 || u[i][1] == 1)
					print(t[0], u[i][3], u[i][2]);
		} else if (indel_bed) {
			var indel_len = [];
			t[1] = parseInt(t[1]);
			for (var i = 0; i < u.length; ++i)
				if (u[i][1] == 3 || u[i][1] == 4) indel_len.push(u[i][4]);
			if (indel_len.length)
				print(t[0], t[1] - 1, t[1] - 1 + t[3].length, indel_len.sort().join(","), u[0][2]);
			continue;
		} else {
			if (u[0][2] < 0) continue;
			for (var i = 0; i < u.length; ++i)
				a.push(u[i].slice(0, 3));
		}
	}
	buf.destroy(buf);
	file.close();

	a.sort(function(x,y) {return y[0]-x[0]});
	if (!indel_bed && !snp_pos) {
		if (show_hdr) print("Q", "#", "#tsHom", "#tsHet", "#tvHom", "#tvHet", "#mnpHom", "#mnpHet", "#insHom", "#insHet", "#delHom", "#delHet");
		var size = Math.floor(a.length * b + 1.);
		var lastq = -1;
		var c = [], ac = [];
		for (var j = 0; j < 11; ++j) c[j] = ac[j] = 0;
		for (var i = 0; i <= a.length; ++i) {
			if (i == a.length || (a[i][0] != lastq && c[0] > a.length * b)) {
				print(lastq, ac.join("\t"), c.join("\t"));
				if (i == a.length) break;
				for (var j = 0; j < 11; ++j) c[j] = 0;
			}
			++c[0]; ++ac[0];
			++c[a[i][1] * 2 + a[i][2] + 1];
			++ac[a[i][1] * 2 + a[i][2] + 1];
			lastq = a[i][0];
		}
	}
}

/************************************************
 *** recompute GQ and GT in single-sample VCF ***
 ************************************************/

function b8_upd1gt(args)
{
	var c, prior = 30, quiet = false;
	while ((c = getopt(args, "p:q")) != null)
		if (c == 'p') prior = parseFloat(getopt.arg);
		else if (c == 'q') quiet = true;

	// each value in GL corresponds to which genotype? Assuming diploid
	var max_alleles = 16, het_arr = [], gt_arr = [];
	for (var i = 0; i < max_alleles; ++i) {
		for (var j = 0; j <= i; ++j) {
			het_arr.push(i != j);
			gt_arr.push(j + '/' + i);
		}
	}

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var re = /^(\d+)(\/|\|)(\d+)/;
	var lineno = 0;
	var have_GQ = false;

	while (file.readline(buf) >= 0) {
		++lineno;
		var line = buf.toString();
		if (line.charAt(0) == '#') {
			if (/^##FORMAT=<ID=GQ/.test(line))
				have_GQ = true;
			if (line.charAt(1) != '#' && !have_GQ)
				print('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">');
			print(line);
			continue;
		}
		var t = line.split("\t");
		var s1 = t[8].split(":");
		var s2 = t[9].split(":");
		var i;
		for (i = 0; i < s1.length; ++i)
			if (s1[i] == 'GL' || s1[i] == 'PL') break;
		if (i == s1.length) {
			if (!quiet) warn("WARNING: no GL/PL at line " + lineno + ":\n" + buf.toString());
			print(buf);
			continue;
		}
		var u = s2[i].split(",");
		var n_alleles = t[4].split(",").length + 1;
		if (u.length != n_alleles * (n_alleles + 1) / 2) {
			if (!quiet) warn("WARNING: inconsistent GL/PL at line " + lineno + ":\n" + buf.toString());
			print(buf);
			continue;
		}
		var is_GL = (s1[i] == 'GL');
		for (var j = 0; j < u.length; ++j) {
			u[j] = parseFloat(u[j]);
			if (is_GL) u[j] *= -10.;
			u[j] += (het_arr[j]? prior : 0);
		}
		var min = 1e37, min_j = -1, min2 = 1e37;
		for (var j = 0; j < u.length; ++j)
			if (min > u[j]) min2 = min, min = u[j], min_j = j;
			else if (min2 > u[j]) min2 = u[j]; 
		var GT = gt_arr[min_j], GQ = Math.floor(min2 - min + .499);
		if (s1[0] != 'GT') {
			s1.unshift('GT'); s2.unshift(GT);
		} else s2[0] = GT;
		var k;
		for (k = 0; k < s1.length; ++k)
			if (s1[k] == 'GQ') break;
		if (k == s1.length) {
			s1.push('GQ'); s2.push(GQ);
		} else s2[k] = GQ;
		t[8] = s1.join(":");
		t[9] = s2.join(":");
		print(t.join("\t"));
	}
	buf.destroy(); file.close();
}

/******************
 *** De-overlap ***
 ******************/

function b8_deovlp(args)
{
	var c;
	while ((c = getopt(args, "")) != null);

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var c1 = '#'.charCodeAt(0);
	var a = [];
	var n = 0;
	while (file.readline(buf) >= 0) {
		if (buf[0] == c1) {
			print(buf);
			continue;
		}
		var s = buf.toString();
		var t = s.split("\t", 6);
		t[1] = parseInt(t[1]) - 1; t[5] = parseFloat(t[5]);
		if (a.length) {
			var i;
			for (i = 0; i < a.length && (a[i][0] != t[0] || a[i][1] <= t[1]); ++i)
				if (a[i][3]) print(a[i][4]);
				else ++n;
			while (i--) a.shift();
		}
		var to_print = true;
		if (a.length) {
			for (var i = 0; i < a.length; ++i) {
				if (a[i][1] <= t[1] || !a[i][3]) continue;
				if (a[i][2] < t[5]) a[i][3] = false;
				else to_print = false;
			}
		}
		a.push([t[0], t[1] + t[3].length, t[5], to_print, s]);
	}
	for (var i = 0; i < a.length; ++i)
		if (a[i][3]) print(a[i][4]);
		else ++n;
	buf.destroy();
	file.close();
	warn(n + " variants have been dropped");
}

/*************************************
 * Convert CG's masterVarBeta to VCF *
 ************************************/

function b8_cg2vcf(args)
{
	var file = args.length > 0? new File(args[0]) : new File();
	var buf = new Bytes();
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var a = [];
		if (t.length < 5 || buf[0] < 48 || buf[0] > 57 || t[5] == 'no-call' || t[6] == 'ref') continue;
		if (t[7].length != t[8].length || t[7].length != t[9].length)
			for (var i = 7; i <= 9; ++i)
				if (t[i] != '?') t[i] = 'N' + t[i];
		var dp = t[24] == ''? 0 : parseInt(t[24]);
		var call = [], dr = 0, da = 0, gq = 1000000, qual = 0;
		for (var i = 0; i < 2; ++i) {
			var q = parseInt(t[10+i]);
			gq = gq < q? gq : q;
			if (t[i+8] == '?') continue;
			call.push(t[i+8]);
			var d = t[21+i] == ''? 0 : parseInt(t[21+i]);
			if (t[i+8] != t[7]) {
				da += i == 0 || t[9] != t[8]? d : 0;
				qual = qual > q? qual : q;
			} else dr += d;
		}
		var gt = './.', alt = '';
		if (call.length == 2) {
			if (call[0] == call[1]) {
				gt = call[0] == t[7]? '0/0' : '1/1';
				alt = call[0] == t[7]? '.' : call[0];
			} else if (call[0] == t[7] || call[1] == t[7]) { // biallelic
				gt = '0/1';
				alt = call[0] != t[7]? call[0] : call[1];
			} else { // triallelic
				gt = '1/2';
				alt = call.join(",");
			}
		} else alt = call[0];
		var o = [t[2], parseInt(t[3]) + 1, '.', t[7], alt, qual, '.', "DP="+dp+";"+"CGDP2="+dr+","+da, "GT:GQ", gt+":"+gq];
		print(o.join("\t"));
	}
	buf.destroy();
	file.close();
}

/********************
 *** Annotate VCF ***
 ********************/

function b8_anno(args)
{
	var c, bed = [], depth_only = false, no_header = false, clear_flt = false, hard_Q = null;
	while ((c = getopt(args, 'db:Q:HFh')) != null) {
		if (c == 'b') bed.push(read_bed(getopt.arg)[0]);
		else if (c == 'd') depth_only = true;
		else if (c == 'Q') hard_Q = parseInt(getopt.arg);
		else if (c == 'H') no_header = true;
		else if (c == 'F') clear_flt = true;
		else if (c == 'h') {
			print("\nUsage:   k8 hapdip.js anno [options] <in.vcf>\n");
			print("Options: -Q INT     drop variants with quality below INT [null]");
			print("         -b FILE    mark variants overlapping BED FILE [null]");
			print("         -H         suppress the VCF header");
			print("         -F         clear the FILTER field in VCF");
			print("");
			exit(1);
		}
	}

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();

	var csharp = '#'.charCodeAt(0);
	var lineno = 0, sum_dp = 0, n_dp = 0;
	while (file.readline(buf) >= 0) {
		++lineno;
		if (buf[0] == csharp) {
			if (!no_header) {
				if (buf[1] != csharp) { // print extra INFO tags
					print('##INFO=<ID=_DP,Number=1,Type=Integer,Description="Raw read depth">');
					print('##INFO=<ID=_DS,Number=1,Type=Integer,Description="min{alt_DP_on_forward, alt_DP_on_reverse}">');
					print('##INFO=<ID=_AB,Number=1,Type=Integer,Description="Percentage of non-reference reads">');
					print('##INFO=<ID=_FS,Number=1,Type=Integer,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">');
					for (var i = 0; i < bed.length; ++i)
						print('##INFO=<ID=_BED'+(i+1)+',Number=1,Type=Integer,Description="Overlapping with BED">');
				}
				print(buf);
			}
			continue;
		}
		var t = buf.toString().split("\t");
		if (hard_Q != null && parseFloat(t[5]) < hard_Q) continue;
		t[1] = parseInt(t[1]);
		// extract depth information
		var m = /DP=(\d+)/.exec(t[7]);
		var depth = m != null? parseInt(m[1]) : -1; // get read depth
		var dp4 = [], dp_ref = null, dp_alt = null, dp_alt_for = null, dp_alt_rev = null, FS = null;
		var m4 = [];
		// get generic depth information
		if (t[9] != null && t[8] != null) {
			var AD = null, ADF = null, ADR = null, DP = null;
			var fmt = t[8].split(":");
			for (var i = 0; i < fmt.length; ++i) {
				if (fmt[i] == 'AD') AD = i;
				else if (fmt[i] == 'ADF') ADF = i;
				else if (fmt[i] == 'ADR') ADR = i;
				else if (fmt[i] == 'DP') DP = i;
			}
			if (ADF != null && ADR != null) {
				var s = t[9].split(":");
				var cf = s[ADF].split(",");
				var cr = s[ADR].split(",");
				dp4 = [parseInt(cf[0]),parseInt(cr[0]),0,0];
				for (var i = 1; i < cf.length; ++i)
					dp4[2] += parseInt(cf[i]), dp4[3] += parseInt(cr[i]);
			} else if (AD != null) {
				var a = t[9].split(":")[AD].split(",");
				dp_ref = parseInt(a[0]), dp_alt = 0;
				for (var i = 1; i < a.length; ++i)
					dp_alt += parseInt(a[i]);
			} else if (DP != null) {
				depth = parseInt(t[9].split(":")[DP]);
			}
		}
		if ((m = /[;\t]FS=([^\t;]+)/.exec(t[7])) != null)
			FS = Math.pow(10, -.1 * parseFloat(m[1]));
		// caller specific depth
		if ((m = /DP4=(\d+),(\d+),(\d+),(\d+)/.exec(t[7])) != null) { // samtools
			for (var j = 1; j <= 4; ++j) dp4[j-1] = parseInt(m[j]);
		} else if ((m = /CGDP2=(\d+),(\d+)/.exec(t[7])) != null) { // CG vcf converted by hapdip.js cg2vcf
			dp_ref = parseInt(m[1]);
			dp_alt = parseInt(m[2]);
		} else if ((m4[0] = /SRF=(\d+)/.exec(t[7])) != null && (m4[1] = /SRR=(\d+)/.exec(t[7])) != null
				&& (m4[2] = /SAF=(\d+)/.exec(t[7])) != null && (m4[3] = /SAR=(\d+)/.exec(t[7])) != null) { // freebayes; in four separate tags
			for (var j = 0; j < 4; ++j) dp4[j] = parseInt(m4[j][1]);
		} else if ((m = /NF=([\d,]+).*NR=([\d,]+).*TCF=(\d+).*TCR=(\d+)/.exec(t[7])) != null) { // Platypus; in four tags, but I don't really know what they mean...
			var m1 = m[1].split(",");
			var m2 = m[2].split(",");
			var max = -1, max_j = -1;
			if (m1.length != m2.length)
				throw Error("Platypus: INFO:NF and INFO:NR are inconsistent");
			for (var j = 0; j < m1.length; ++j) {
				m1[j] = parseInt(m1[j]); m2[j] = parseInt(m2[j]);
				var x = m1[j] + m2[j]
					if (max < x) max = x, max_j = j;
			}
			dp4[2] = m1[max_j]; dp4[3] = m2[max_j];
			dp4[0] = parseInt(m[3]) - dp4[2];
			dp4[1] = parseInt(m[4]) - dp4[3];
			if (dp4[0] < 0 || dp4[1] < 0) {
				warn("Platypus: negative reference coverage at line " + lineno);
				if (dp4[0] < 0) dp4[0] = 0;
				if (dp4[1] < 0) dp4[1] = 0;
			}
		}

		if (dp4.length) dp_ref = dp4[0] + dp4[1], dp_alt_for = dp4[2], dp_alt_rev = dp4[3], dp_alt = dp_alt_for + dp_alt_rev;
		if (dp_alt != null && dp_ref != null) depth = dp_alt + dp_ref;
		if (depth == null) throw Error("Cannot get the depth information at line " + lineno + ":\n" + buf.toString());
		if (FS == null && dp4.length) FS = Math.fisher_exact(dp4[0], dp4[1], dp4[2], dp4[3])[0];
		sum_dp += depth; ++n_dp;

		var labels = ['_DP', '_DS', '_AB', '_FS'];
		//             0      1      2      3
		var values = [null,   null,  null,  null];
		values[0] = depth; // DP: total depth
		if (dp_alt_for != null && dp_alt_rev != null) // DS: double-strand support
			values[1] = dp_alt_for < dp_alt_rev? dp_alt_for : dp_alt_rev;
		if (dp_alt != null) // AB: allele balance
			values[2] = depth == 0? 0 : Math.floor(100 * dp_alt / depth + .499);
		if (FS != null) // FS: fisher strand bias
			values[3] = Math.floor((FS < 1e-10? 0 : -4.343 * Math.log(FS)) + .499);
		for (var i = 0; i < bed.length; ++i) {
			var bi = bed[i];
			labels.push('_BED' + (i+1));
			if (bi[t[0]] != null) {
				var x = bed[i][t[0]](t[1] - 1, t[1] - 1 + t[3].length);
				values[4+i] = x == null? 0 : 1;
			}
		}
		// write INFO
		var extra_info = [];
		for (var i = 0; i < values.length; ++i)
			if (values[i] != null)
				extra_info.push(labels[i] + '=' + values[i]);
		t[7] = extra_info.join(';') + (t[7] == "." ? '' : ';' + t[7]);
		if (clear_flt) t[6] = '.';
		if (depth_only) print(t[0], t[1], depth);
		else print(t.join("\t"));
	}

	buf.destroy();
	file.close();
	warn("Average depth: " + (sum_dp/n_dp).toFixed(3));
	warn("Suggested depth cut-off: " + (sum_dp/n_dp + 4*Math.sqrt(sum_dp/n_dp)).toFixed(3));
}

/*************************
 *** Filter anno'd VCF ***
 *************************/

function b8_filter(args)
{
	var c, AB = 30, AB_gap = null, LC = 0, DP_coef = 4, DS = 1, FS = 30, min_q = 30, min_dp = 3, no_header = false, drop_flt = false, auto_only = false, meanDP = -1;
	while ((c = getopt(args, "Aa:b:q:f:HDd:c:s:F:")) != null) {
		if (c == 'a') AB = parseInt(getopt.arg);
		else if (c == 'b') AB_gap = parseInt(getopt.arg);
		else if (c == 'q') min_q = parseFloat(getopt.arg);
		else if (c == 'f') FS = parseFloat(getopt.arg);
		else if (c == 'H') no_header = true;
		else if (c == 'D') drop_flt = true;
		else if (c == 'd') min_dp = parseInt(getopt.arg);
		else if (c == 'c') DP_coef = parseFloat(getopt.arg);
		else if (c == 's') DS = parseInt(getopt.arg);
		else if (c == 'A') auto_only = true;
		else if (c == 'F') meanDP = parseFloat(getopt.arg);
	}
	if (AB_gap == null) AB_gap = AB;

	if (getopt.ind + 1 > args.length) {
		print("\nUsage:   k8 hapdip.js filter [options] <anno.vcf>\n");
		print("Options: -a INT     min _AB at SNPs ["+AB+"]");
		print("         -b INT     min _AB at INDELs [same as -a]");
		print("         -q FLOAT   min QUAL ["+min_q+"]");
		print("         -f FLOAT   max _FS ["+FS+"]");
		print("         -d INT     min _DP ["+min_dp+"]");
		print("         -F FLOAT   set meanDP to FLOAT [inferred]");
		print("         -c FLOAT   set max _DP to 'meanDP + FLOAT * sqrt(meanDP)' ["+DP_coef+"]");
		print("         -s INT     min _DS ["+DS+"]");
		print("         -H         suppress the VCF header");
		print("         -D         drop filtered variants");
		print("         -A         only output lines matching /^(chr)?[0-9]+/");
		print("");
		exit(1);
	}

	var file = new File(args[getopt.ind]);
	var buf = new Bytes();

	// pass 1: compute avg depth
	var sum_dp = 0, n_dp = 0;
	warn("Pass 1: computing average depth at variant sites...");
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (auto_only && !/^(chr)?[0-9]+$/.test(t[0])) continue;
		var m;
		if (parseFloat(t[5]) < min_q) continue;
		if ((m = /_DP=(\d+)/.exec(t[7])) == null) continue;
		var dp = parseInt(m[1]);
		sum_dp += dp; ++n_dp;
	}
	var avg_dp = sum_dp / n_dp;
	if (meanDP > 0) avg_dp = meanDP;
	var max_dp = Math.sqrt(avg_dp) * DP_coef + avg_dp;
	var DP = max_dp;
	file.close();
	warn("Average depth: " + avg_dp.toFixed(2));
	warn("Max depth cutoff: " + max_dp.toFixed(2));

	// pass 2: filter
	warn("Pass 2: applying filters...");
	file = new File(args[getopt.ind]);
	var max_BED = 0;
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (/^#/.test(line)) {
			var m;
			if ((m = /INFO.*ID=_BED(\d+)/.exec(line)) != null)
				max_BED = max_BED > parseInt(m[1])? max_BED : parseInt(m[1]);
			if (no_header) continue;
			if (/^#[^#]/.test(line)) {
				print('##FILTER=<ID=DPhigh,Description="High read depth: _DP>'+max_dp.toFixed(2)+'">');
				print('##FILTER=<ID=DPlow,Description="Low read depth: _DP<'+min_dp+'">');
				print('##FILTER=<ID=FShigh,Description="Large Fisher-Strand bias: _FS>'+FS+'">');
				print('##FILTER=<ID=ABlow,Description="Low fraction of non-reference reads: _AB<'+AB+' at SNPs or _AB<'+AB_gap+' at INDELs">');
				print('##FILTER=<ID=DSlow,Description="Low double-strand support at SNPs: _DS<'+DS+'">');
				for (var i = 0; i < max_BED; ++i)
					print('##FILTER=<ID=BED'+(i+1)+',Description="Overlapping _BED'+(i+1)+'">');
			}
			print(line);
			continue;
		}
		var t = line.split("\t");
		if (auto_only && !/^(chr)?[0-9]+$/.test(t[0])) continue;
		var m, flt = '', dp = null;
		if (parseFloat(t[5]) < min_q) continue;
		if ((m = /_DP=(-?\d+)/.exec(t[7])) != null) dp = parseInt(m[1]);
		if (dp != null) {
			if (dp > DP) flt += "DPhigh;";
			else if (dp < min_dp) flt += "DPlow;";
		} else flt += "noDP";
		if ((m = /_FS=(\d+)/.exec(t[7])) != null && parseInt(m[1]) > FS) flt += "FShigh;";
		var has_gap = false, s = t[4].split(",");
		for (var i = 0; i < s.length; ++i)
			if (s[i].length != t[3].length) has_gap = true;
		if (!has_gap) {
			if ((m = /_DS=(\d+)/.exec(t[7])) != null && parseInt(m[1]) < DS) flt += "DSlow;";
			if ((m = /_AB=(\d+)/.exec(t[7])) != null && parseInt(m[1]) < AB) flt += "ABlow;";
		} else {
			if ((m = /_AB=(\d+)/.exec(t[7])) != null && parseInt(m[1]) < AB_gap) flt += "ABlow;";
		}
		if ((m = /_BED(\d+)=(-?\d+)/.exec(t[7])) != null && parseInt(m[2]) != 0) flt += "BED"+m[1]+';';
		if (drop_flt && flt != '') continue;
		t[6] = flt == ''? "." : flt.substr(0, flt.length-1);
		print(t.join("\t"));
	}
	file.close();

	buf.destroy();
}

/********************************
 *** Convert VCF to unary BED ***
 ********************************/

// the following is similar to b8_parse_vcf1
function b8_parse_vcf2(t) // t = vcf_line.split("\t")
{
	if (t.length < 6) return null;
	var a = [], pos = parseInt(t[1]) - 1;
	t[3] = t[3].toUpperCase(); t[4] = t[4].toUpperCase();
	var s = t[4].split(","); // list of ALT alleles
	// get CIGAR for freebayes
	var m3 = /CIGAR=([^;\t]+)/.exec(t[7]);
	var cigar = m3 != null? m3[1].split(",") : [];
	if (cigar.length && cigar.length != s.length) throw Error("Inconsistent ALT and CIGAR");
	// loop through each ALT allele
	for (var i = 0; i < s.length; ++i) {
		if (t[3].length == 1 && s[i].length == 1) { // SNP
			if (t[3] != s[i]) a.push([pos, pos+1, 0, t[3], s[i], i]);
		} else if (cigar.length) { // MNP or INDEL from freebayes
			var x = 0, y = 0;
			var m4, re = /(\d+)([MXID])/g;
			while ((m4 = re.exec(cigar[i])) != null) {
				var l = parseInt(m4[1]);
				if (m4[2] == 'X') {
					for (var j = 0; j < l; ++j) {
						var u = t[3].charAt(x+j), v = s[i].charAt(y+j);
						a.push([pos + x, pos+x+1, 0, u, v, i]);
					}
					x += l, y += l;
				} else if (m4[2] == 'I') {
					if (x == 0 || y == 0) throw Error("Leading I/D");
					var u = t[3].substr(x-1, 1), v = s[i].substr(y-1, l+1);
					a.push([pos + x - 1, pos + x, l, u, v, i]);
					y += l;
				} else if (m4[2] == 'D') {
					if (x == 0 || y == 0) throw Error("Leading I/D");
					var u = t[3].substr(x-1, l+1), v = s[i].substr(y-1, 1);
					a.push([pos + x - 1, pos + x + l, -l, u, v, i]);
					x += l;
				} else if (m4[2] == 'M') x += l, y += l;
			}
		} else { // MNP or INDEL from Platypus and others; FIXME: use the b8_parse_vcf1() approach here!
			var st = b8_allele_start(t[3], s[i]);
			var d = s[i].length - t[3].length, rs, as;
			if (d > 0) {
				a.push([pos + st - 1, pos + st, d, t[3].charAt(st - 1), s[i].substr(st - 1, d + 1), i]);
				rs = st, as = d + st;
			} else if (d < 0) {
				a.push([pos + st - 1, pos + st + (t[3].length + d), d, t[3].substr(st - 1, -d + 1), s[i].charAt(st - 1), i]);
				rs = -d + st, as = st;
			} else rs = as = 0;
			for (var j = 0; j < t[3].length - rs; ++j) { // check the rest
				var u = t[3].charAt(rs + j), v = s[i].charAt(as + j);
				if (u != v) a.push([pos + j + rs, pos + j + rs + 1, 0, u, v, i]);
			}
		}
	}
	return a; // [start, end, indelLen, ref, alt, i]
}

function b8_get_ACA(t)
{
	var n_alleles = t[4].split(",").length + 1;
	var ACA = [];
	for (var i = 0; i < n_alleles; ++i) ACA[i] = 0;
	if (t.length >= 10 && t[8].substr(0, 2) == 'GT') {
		for (var i = 9; i < t.length; ++i) {
			var m;
			if ((m = /(\.|\d+)[\/\|](\.|\d+)/.exec(t[i])) != null) {
				if (m[1] != '.') ++ACA[parseInt(m[1])];
				if (m[2] != '.') ++ACA[parseInt(m[2])];
			}
		}
	}
	return ACA;
}

function b8_vcf2bed(args)
{
	var file = args.length == 0? new File() : new File(args[0]);
	var buf = new Bytes();

	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '#') continue;
		var t = line.split("\t");
		var a = b8_parse_vcf2(t);
		var ACA = b8_get_ACA(t), tot = 0;
		for (var i = 0; i < ACA.length; ++i) tot += ACA[i];
		for (var i = 0; i < a.length; ++i) {
			a[i][5] = ACA[a[i][5]+1]; a[i][6] = tot;
			print(t[0], a[i].join("\t"), t[6]);
		}
	}

	buf.destroy();
	file.close();
}

/*******************************************
 *** subset, replace and reorder samples ***
 *******************************************/

function b8_vcfsub(args)
{
	if (args.length < 2) {
		print("Usage: k8 vcfsub <list.txt> <in.vcf>");
		exit(1);
	}
	var buf = new Bytes();
	var file = new File(args[0]);
	var list_ori = [];
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t.length >= 2) list_ori.push([t[0], t[1]]);
		else list_ori.push([t[0], null]);
	}
	file.close();

	file = new File(args[1]);
	var map = [], list = [], replace = true;
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '#') {
			if (line.charAt(1) != '#') {
				var t = line.split("\t"), hash = {};
				// get the intersection of list_ori[] and the list in VCF
				for (var i = 9; i < t.length; ++i)
					hash[t[i]] = i;
				for (var j = 0; j < list_ori.length; ++j)
					if (hash[list_ori[j][0]] != null)
						list.push(list_ori[j]);
				// generate map[]
				for (var j = 0; j < list.length; ++j)
					map.push(hash[list[j][0]]);
				// reorder samples
				var s = t.slice(0, 9);
				for (var j = 0; j < list.length; ++j)
					s.push(list[j][list[j][1] == null? 0 : 1]);
				print(s.join("\t"));
			} else print(line);
		} else {
			var t = line.split("\t"), s = t.slice(0, 9);
			for (var j = 0; j < map.length; ++j)
				s.push(t[map[j]]);
			print(s.join("\t"));
		}
	}
	file.close();
	buf.destroy();
}

function b8_vcfsum(args)
{
	var c, beds = [], min_top_sr = 10, no_filter = true, no_gt = false, max_no_call = 0.5;
	while ((c = getopt(args, "Gfb:n:c:")) != null) {
		if (c == 'n') min_top_sr = parseInt(getopt.arg);
		else if (c == 'f') no_filter = false;
		else if (c == 'G') no_gt = true;
		else if (c == 'c') max_no_call = parseFloat(getopt.arg);
		else if (c == 'b') {
			var m, label = null, fn = null;
			if ((m = /([^\s=]+)=(\S+)/.exec(getopt.arg)) != null) {
				label = m[1];
				fn = m[2];
			} else label = fn = getopt.arg;
			warn('Reading BED file ' + fn + '...');
			beds.push([label, read_bed(fn, 11)[0]]);
		}
	}

	if (args.length == getopt.ind) {
		print("Usage: k8 hapdip.js vcfsum [options] <htsbox-pileup.vcf>");
		print("Options:");
		print("  -b STR=FILE   flag STR in INFO if overlapping regions in BED FILE [null]");
		print("  -n INT        threshold on the SR filter (effective with -f) ["+min_top_sr+"]");
		print("  -c FLOAT      threshold on the no-call rate ["+max_no_call+"]");
		print("  -G            discard genotype fields");
		print("  -f            add HIGHNC, LOWAD and NONVAR to FILTER");
		exit(1);
	}

	warn('Processing VCF...');

	var new_lines = "";
	new_lines += '##FILTER=<ID=LOWAD,Description="TOPAD < ' + min_top_sr + '">\n';
	new_lines += '##FILTER=<ID=HIGHNC,Description="no-call rate > ' + max_no_call +'">\n';
	new_lines += '##FILTER=<ID=NONVAR,Description="not a variant site">\n';
	for (var i = 0; i < beds.length; ++i)
		new_lines += '##INFO=<ID=' + beds[i][0] + ',Number=0,Type=Flag>\n';
	new_lines += '##INFO=<ID=TOPAD,Number=R,Type=Integer,Description="best AD across samples">\n';
	new_lines += '##INFO=<ID=AD,Number=R,Type=Integer,Description="allelic depths">\n';
	new_lines += '##INFO=<ID=CA,Number=R,Type=Integer,Description="count of reference and alternate alleles">\n';
	new_lines += '##INFO=<ID=CG,Number=G,Type=Integer,Description="count of genotypes">';

	var file = new File(args[getopt.ind]);
	var buf = new Bytes();

	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '#') {
			if (line.charAt(1) != '#') {
				print(new_lines);
				if (no_gt) line = line.split("\t", 8).join("\t");
			}
			print(line);
			continue;
		}
		var t = line.split("\t");
		var gt = null, sr = null;
		var s = t[8].split(":");
		t[7] = t[7].replace(/\b(TOPAD|AD|CA|CG)=[^\s;]+(;?)/g, '');
		for (var i = 0; i < s.length; ++i)
			if (s[i] == 'GT') gt = i;
			else if (s[i] == 'AD') sr = i;
		var ACA = [], SR = [], topSR = [], CG = []; // for historical reasons, "SR"=="AD"
		var n_alleles = t[4].split(",").length + 1, n_gt = Math.floor(n_alleles * (n_alleles + 1) / 2 + .499);
		for (var i = 0; i < n_alleles; ++i) ACA[i] = SR[i] = topSR[i] = 0;
		for (var i = 0; i < n_gt; ++i) CG[i] = 0;
		for (var j = 9; j < t.length; ++j) {
			s = t[j].split(":");
			if (gt != null) {
				var m;
				if ((m = /(\.|\d+)[\/\|](\.|\d+)/.exec(s[gt])) != null) {
					var h1 = m[1] != '.'? parseInt(m[1]) : -1;
					var h2 = m[2] != '.'? parseInt(m[2]) : -1;
					if (h1 >= 0) ++ACA[h1];
					if (h2 >= 0) ++ACA[h2];
					if (h1 >= 0 && h2 >= 0) {
						var tmp;
						if (h1 > h2) tmp = h1, h1 = h2, h2 = tmp;
						tmp = Math.floor(h1 * (h1 + 1) / 2 + h2 + .499);
						++CG[tmp];
					}
				}
			}
			if (sr != null) {
				var v = s[sr].split(",");
				for (var k = 0; k < v.length; ++k) {
					var u = parseInt(v[k]);
					SR[k] += u;
					topSR[k] = topSR[k] > u? topSR[k] : u;
				}
			}
		}
		// set INFO
		if (t[7] == '.') t[7] = '';
		else if (t[7] != '') t[7] += ';';
		t[7] += 'CA=' + ACA.join(",") + ';CG=' + CG.join(",");
		if (sr != null) t[7] += ';AD=' + SR.join(",") + ';TOPAD=' + topSR.join(",");
		// test BED
		for (var i = 0; i < beds.length; ++i) {
			var start = parseInt(t[1]) - 1;
			if (beds[i][1][t[0]] && beds[i][1][t[0]](start, start + t[3].length) != null)
				t[7] += ';' + beds[i][0];
		}
		// set FILTER
		if (!no_filter) {
			var alt_cnt = 0, flt = '', max_topSR = 0;
			for (var i = 1; i < ACA.length; ++i)
				alt_cnt += ACA[i], max_topSR = max_topSR > topSR[i]? max_topSR : topSR[i];
			if (alt_cnt == 0) flt = 'NONVAR';
			if (sr != null && max_topSR < min_top_sr)
				flt = flt == ''? 'LOWAD' : flt + ';LOWAD';
			if (ACA[0] + alt_cnt < 2 * (t.length - 9) * max_no_call)
				flt = flt == ''? 'HIGHNC' : flt + ';HIGHNC';
			if (t[6] == '.' || t[6] == 'PASS') t[6] = flt == ''? 'PASS' : flt;
			else if (flt != '') t[6] += ';' + flt;
		}
		// print out
		if (no_gt) t.length = 8;
		print(t.join("\t"));
	}

	buf.destroy();
	file.close();
}

function b8_vcfswap(args)
{
	var file = args.length == 0? new File() : new File(args[0]);
	var buf = new Bytes();

	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '#') {
			print(line);
			continue;
		}
		var m, t = line.split("\t");
		for (var i = 9; i < t.length; ++i) {
			if ((m = /^(\d+)([\|\/])(\d+)(.*)/.exec(t[i])) != null) {
				if (m[2] == '|') continue;
				var a = parseInt(m[1]), b = parseInt(m[3]);
				if (a <= b) continue;
				if (m[4]) t[i] = b + "/" + a + m[4];
				else t[i] = b + "/" + a;
			}
		}
		print(t.join("\t"));
	}

	buf.destroy();
	file.close();
}

/**********************************
 *** Evaluate CHM1-NA12878 VCFs ***
 **********************************/

function b8_eval(args)
{
	var c, label = "hapdip", in_bed = [], out_bed = [], drop_flt = true, auto_only = false, min_q = 0;
	while ((c = getopt(args, "aFb:B:L:q:")) != null) {
		if (c == 'b') in_bed.push(read_bed(getopt.arg)[0]);
		else if (c == 'q') min_q = parseInt(getopt.arg);
		else if (c == 'L') label = getopt.arg;
		else if (c == 'B') out_bed.push(read_bed(getopt.arg)[0]);
		else if (c == 'F') drop_flt = false;
		else if (c == 'a') auto_only = true;
	}
	if (getopt.ind + 2 > args.length) {
		print("\nUsage:   k8 hapdip.js [options] <CHM1.vcf> <NA12878.vcf>\n");
		print("Options: -F        ignore the FILTER field in VCF");
		print("         -a        only evaluate lines matching /^(chr)?[0-9]+/");
		print("         -q INT    min base quality");
		print("         -b FILE   only evaluate variants overlapping BED FILE (can be multi) [null]");
		print("         -B FILE   exclude variants overlapping BED FILE (can be multi) [null]");
		print("         -L STR    change label in the output to STR [hapdip]");
		print("");
		exit(1);
	}

	function read_vcf(fn)
	{
		warn("Processing file "+fn+"...");
		var file = new File(fn);
		var buf = new Bytes();
		var cnt = [0, 0, 0]; // SNP, 1bp indel, >1bp indel
		var sharp = '#'.charCodeAt(0);

		while (file.readline(buf) >= 0) {
			if (buf[0] == sharp) continue; // skip header lines
			var t = buf.toString().split("\t");
			if (auto_only && !/^(chr)?[0-9]+$/.test(t[0])) continue;
			if (parseInt(t[5]) < min_q) continue;
			if (drop_flt && (t[6] != "." && t[6] != "PASS")) continue;
			if (in_bed.length || out_bed.length) {
				var start = parseInt(t[1]);
				var end = start + t[3].length;
				var to_proceed = true;
				for (var i = 0; i < in_bed.length; ++i) {
					if (in_bed[i][t[0]] == null || in_bed[i][t[0]](start, end) == null) {
						to_proceed = false;
						break;
					}
				}
				for (var i = 0; i < out_bed.length; ++i) {
					if (out_bed[i][t[0]] != null && out_bed[i][t[0]](start, end) != null) {
						to_proceed = false;
						break;
					}
				}
				if (!to_proceed) continue;
			}
			var c = b8_parse_vcf1(t);
			for (var i = 0; i < c.length; ++i) {
				if (c[i][2] == 0) continue;
				if (c[i][1] == 0 || c[i][1] == 1) ++cnt[0];
				else if (c[i][1] == 3 || c[i][1] == 4) {
					if (c[i][4] == 1 || c[i][4] == -1) ++cnt[1];
					else ++cnt[2];
				}
			}
		}

		buf.destroy();
		file.close();
		return cnt;
	}

	var hetcnt = [];
	hetcnt[0] = read_vcf(args[getopt.ind]);
	hetcnt[1] = read_vcf(args[getopt.ind+1]);
	print(label, "SNP", "FP", hetcnt[0][0]);
	print(label, "SNP", "TP", hetcnt[1][0] - hetcnt[0][0]);
	print(label, "INDEL", "FP", hetcnt[0][1] + hetcnt[0][2]);
	print(label, "INDEL", "TP", (hetcnt[1][1] + hetcnt[1][2]) - (hetcnt[0][1] + hetcnt[0][2]));
	print(label, "INDEL", "FP", hetcnt[0][1], "INDEL-1bp");
	print(label, "INDEL", "TP", hetcnt[1][1] - hetcnt[0][1], "INDEL-1bp");
}

/*****************************
 * Distance based evaluation *
 *****************************/

function b8_read_hrun(fn_hp)
{
	var file = new File(fn_hp);
	var buf = new Bytes();
	var hp = new Map();
	warn("Reading homopolymer runs...");
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t.length < 4) continue;
		hp.put(t[0] + ':' + t[1] + t[3].toUpperCase(), t[2] - t[1]);
	}
	buf.destroy();
	file.close();
	return hp;
}

function b8_merge_print_mt(s, t)
{
	if (s[1] == t[1] && s[3].length > t[3].length) {
		var x;
		x = s, s = t, t = x;
	}
	var m, l = t[1] - s[1], r = l + t[3].length - s[3].length;
	var left = s[3].substr(0, l);
	var right = t[3].substr(t[3].length - r);
	var ref = s[3] + right;
	var alt = s[4] + right + "," + left + t[4];
	s[3] = ref, s[4] = alt, s[9] = "1/2";
	var hpl_s = (m = /HPL=(\d+)/.exec(s[7])) != null? m[1] : '.';
	var hpl_t = (m = /HPL=(\d+)/.exec(t[7])) != null? m[1] : '.';
	if (hpl_s != '.' || hpl_t != '.')
		s[7] = "HPL=" + hpl_s + "," + hpl_t;
	print(s.join("\t"));
}

function b8_postatom(args)
{
	var fn_hp = null;
	while ((c = getopt(args, "p:")) != null)
		if (c == 'p') fn_hp = getopt.arg;
	if (getopt.ind + 1 > args.length) {
		print("Usage: bgt atomize -S <in.vcf> | k8 hapdip.js postatom [-p hrun] /dev/stdin");
		exit(1);
	}

	var hp = fn_hp != null? b8_read_hrun(fn_hp) : null;
	var buf = new Bytes();
	warn("Processing VCF...");
	var file = new File(args[getopt.ind]);
	var last_mt = null, last_hrun = 0;
	var ctg_dict = {};
	while (file.readline(buf) >= 0) {
		var m, l = buf.toString();
		if (l.charAt(0) == '#') {
			if ((m = /^##contig.*ID=([^,\s]+)/.exec(l)) != null) {
				if (ctg_dict[m[1]] != null) // skip duplicated contig name
					continue;
				ctg_dict[m[1]] = 1;
			}
			if (hp != null && l.charAt(1) != '#')
				print('##INFO=<ID=HPL,Number=A,Type=Integer,Description="Homopolymer run length, if available">');
			print(l);
		} else {
			var t = l.split("\t");
			t[1] = parseInt(t[1]);
			var pos = t[1] - 1;
			var is_mt = (/^\.(\/|\|)1/.test(t[9]) || /^1(\/|\|)\./.test(t[9]));
			var gap = t[4].length - t[3].length;
			var hrun = 0;
			if (hp != null && gap != 0) { // test if we are at a homopolymer INDEL
				var j, seq = gap > 0? t[4].substr(1) : t[3].substr(1);
				for (j = 1; j < seq.length; ++j) // test if ins/del bases are all the same
					if (seq.charAt(j) != seq.charAt(0))
						break;
				if (j == seq.length) {
					var key = t[0] + ':' + (pos+1) + seq.charAt(0);
					var val = hp.get(key);
					if (val != null) // a homopolymer indel
						hrun = parseInt(val);
				}
			}
			if (hp != null && hrun > 0) t[7] = "HPL=" + hrun;
			if (last_mt != null) {
				if (!is_mt || (t[0] == last_mt[0] && pos >= last_mt[1]-1 + last_mt[3].length)) {
					print(last_mt.join("\t"));
					last_mt = null;
				} else {
					b8_merge_print_mt(last_mt, t);
					last_mt = null;
					continue;
				}
			}
			if (is_mt) last_mt = t.slice(0);
			else if (hrun > 0) print(t.join("\t"));
			else print(l);
		}
	}
	if (last_mt != null) print(last_mt.join("\t"));
	file.close();
	buf.destroy();
	if (hp != null) hp.destroy();
}

function b8_atomcnt(args)
{
	var c, min_l = 0, max_l = 1<<30, max_hpl = 1<<30, fn_bed = null, bed = null, to_print = false;
	while ((c = getopt(args, "b:l:L:p:P")) != null) {
		if (c == 'b') fn_bed = getopt.arg;
		else if (c == 'l') min_l = parseInt(getopt.arg);
		else if (c == 'L') max_l = parseInt(getopt.arg);
		else if (c == 'p') max_hpl = parseInt(getopt.arg);
		else if (c == 'P') to_print = true;
	}
	if (getopt.ind + 1 > args.length) {
		print("Usage: k8 atomcnt [options] <in.vcf>");
		print("Options:");
		print("  -b FILE     ignore variants not overlapping BED FILE [null]");
		print("  -l INT      ignore INDELs shorter than INT [0]");
		print("  -L INT      ignore INDELs longer than INT [inf]");
		print("  -p INT      ignore HPL longer than INT [inf]");
		print("  -P          print counted VCF (don't print counts)");
		exit(1);
	}
	if (fn_bed != null) bed = read_bed(fn_bed)[0];

	var file = new File(args[getopt.ind]);
	var buf = new Bytes();
	var n_sub = 0, n_gap = 0;
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '#') {
			if (to_print) print(line);
			continue;
		}
		var t = line.split("\t");
		var pos = parseInt(t[1]) - 1;
		if (bed != null) {
			if (bed[t[0]] == null) continue;
			if (bed[t[0]](pos, pos + t[3].length) == null) continue;
		}
		var s = t[4].split(",");
		var m, hpl = null, e_sub = false, e_gap = false;
		if ((m = /\bHPL=([^\s;]+)/.exec(t[7])) != null)
			hpl = m[1].split(",");
		for (var i = 0; i < s.length; ++i) {
			var gap = s[i].length - t[3].length;
			if (gap != 0) {
				var gap_len = gap > 0? gap : -gap;
				if (gap_len < min_l || gap_len > max_l) continue;
				if (hpl != null && hpl[i] != null && hpl[i] != '.') {
					var l = parseInt(hpl[i]);
					if (l > max_hpl) continue;
				}
				e_gap = true;
			} else e_sub = true;
		}
		if (e_gap) ++n_gap;
		if (e_sub) ++n_sub;
		if (to_print && (e_gap || e_sub))
			print(line);
	}
	buf.destroy();
	file.close();
	if (!to_print) print(n_sub, n_gap);
}

function b8_distEval(args)
{
	var c, fn_bed = null, fn_hp = null, max_d = 10, min_l = 0, max_l = 1<<30, eval_snp = true, eval_indel = true, show_err = false;
	while ((c = getopt(args, "b:d:ISep:l:L:")) != null) {
		if (c == 'b') fn_bed = getopt.arg;
		else if (c == 'p') fn_hp = getopt.arg;
		else if (c == 'd') max_d = parseInt(getopt.arg);
		else if (c == 'S') eval_snp = false;
		else if (c == 'I') eval_indel = false;
		else if (c == 'e') show_err = true;
		else if (c == 'l') min_l = parseInt(getopt.arg);
		else if (c == 'L') max_l = parseInt(getopt.arg);
	}
	if (getopt.ind + 2 > args.length) {
		print("");
		print("Usage:   k8 hapdip.js distEval [options] <P.vcf> <call.vcf>\n");
		print("Options: -d INT     max distance ["+max_d+"]");
		print("         -b FILE    N+P regions, required unless P.vcf is a gVCF [null]");
		print("         -p FILE    homopolymer BED generated by 'seqtk hrun'");
		print("         -l INT     min INDEL length ["+min_l+"]");
		print("         -L INT     max INDEL length [inf]");
		print("         -S         skip SNPs");
		print("         -I         skip INDELs");
		print("         -e         print FN/FP (fmt: chr, start, end, indel, FN/FP)");
		print("");
		print("Note: By default, if a called SNP is close to a true INDEL but no other");
		print("      true SNPs, it is still considered to be correct. When -S or -I is");
		print("      applied, the above case is counted as an error.");
		print("");
		exit(1);
	}

	var sum_NP = null, bed_NP = null;
	if (fn_bed != null) {
		warn("Reading N+P regions...");
		var a = read_bed(fn_bed);
		sum_NP = a[2], bed_NP = a[0];
	}
	var hp = fn_hp != null? b8_read_hrun(fn_hp) : null;

	function read_vcf(fn)
	{
		var file = new File(fn);
		var buf = new Bytes();
		var sharp = '#'.charCodeAt(0);
		var reg = {}, idx = {};

		while (file.readline(buf) >= 0) {
			if (buf[0] == sharp) continue; // skip header lines
			var t = buf.toString().split("\t");
			if (t[6] != '.' && t[6] != 'PASS') continue; // skip filtered variants
			var a = b8_parse_vcf1(t);
			if (reg[t[0]] == null) reg[t[0]] = [];
			for (var i = 0; i < a.length; ++i) {
				var x = a[i], end, flt = false;
				if (!eval_indel && x[4] != 0) continue;
				if (!eval_snp   && x[4] == 0) continue;
				if (x[4] != 0) {
					if (x[4] < min_l && x[4] > -min_l) flt = true;
					if (x[4] > max_l || x[4] < -max_l) flt = true;
				}
				if (hp != null && x[4] != 0) { // test if we are at a homopolymer INDEL
					var j, seq = x[5];
					for (j = 1; j < seq.length; ++j) // test if ins/del bases are all the same
						if (seq.charAt(j) != seq.charAt(0))
							break;
					if (j == seq.length) {
						var key = t[0] + ':' + x[3] + seq.charAt(0);
						if (hp.get(key) != null) flt = true;
					}
				}
				if (x[4] > 0) reg[t[0]].push([x[3] - 1, x[3] + 1, x[4], flt]);
				else if (x[4] == 0) reg[t[0]].push([x[3], x[3] + 1, x[4], flt]);
				else reg[t[0]].push([x[3], x[3] - x[4], x[4], flt]);
			}
		}

		buf.destroy();
		file.close();
		for (var chr in reg)
			if (reg[chr].length > 0)
				idx[chr] = intv_ovlp3(reg[chr], 13, false)[0];
		return [reg, idx];
	}

	warn("Reading "+args[getopt.ind]+"...");
	var truth = read_vcf(args[getopt.ind]);
	warn("Reading "+args[getopt.ind+1]+"...");
	var call  = read_vcf(args[getopt.ind+1]);
	var TP = [0,0], FP = [0,0], FN = [0,0], TP1 = [0,0];

	warn("Counting TP and FN...");
	var cc = [0, 0, 0, 0, 0];
	for (var chr in truth[0]) {
		var chr_NP = bed_NP != null? bed_NP[chr] : null;
		var chr_call = call[1][chr];
		if (bed_NP != null && chr_NP == null) continue; // not in N+P
		var x = truth[0][chr];
		for (var i = 0; i < x.length; ++i) {
			if (bed_NP != null && chr_NP(x[i][0], x[i][1]) == null) continue; // not in N+P
			if (x[i][3]) continue; // filtered
			var start = x[i][0] - max_d, end = x[i][1] + max_d, type = x[i][2] == 0? 0 : 1;
			if (start < 0) start = 0;
			if (chr_call == null || chr_call(start, end) == null) {
				++FN[type];
				if (show_err) print(chr, x[i][0], x[i][1], x[i][2], 'FN');
			} else ++TP[type];
		}
	}

	warn("Counting FP...");
	for (var chr in call[0]) {
		var chr_NP = bed_NP != null? bed_NP[chr] : null;
		var chr_truth = truth[1][chr];
		if (bed_NP != null && chr_NP == null) continue; // not in N+P
		var x = call[0][chr];
		for (var i = 0; i < x.length; ++i) {
			if (bed_NP != null && chr_NP(x[i][0], x[i][1]) == null) continue; // not in N+P
			if (x[i][3]) continue; // filtered
			var start = x[i][0] - max_d, end = x[i][1] + max_d, type = x[i][2] == 0? 0 : 1;
			if (start < 0) start = 0;
			if (chr_truth == null || chr_truth(start, end) == null) {
				++FP[type];
				if (show_err) print(chr, x[i][0], x[i][1], x[i][2], 'FP');
			} else ++TP1[type];
		}
	}

	if (!show_err) {
		if (eval_snp) {
			if (bed_NP != null)
				print("distEval", 'SNP', "N+P", sum_NP);
			print("distEval", 'SNP', "TP", TP[0]);
			print("distEval", 'SNP', "FN", FN[0]);
			print("distEval", 'SNP', "FP", FP[0]);
			print("distEval", 'SNP', '%FNR', (100 * FN[0] / (FN[0] + TP[0])).toFixed(2));
			print("distEval", 'SNP', '%FDR', (100 * FP[0] / (FP[0] + TP1[0])).toFixed(2));
		}
		if (eval_indel) {
			if (bed_NP != null)
				print("distEval", 'INDEL', "N+P", sum_NP);
			print("distEval", 'INDEL', "TP", TP[1]);
			print("distEval", 'INDEL', "FN", FN[1]);
			print("distEval", 'INDEL', "FP", FP[1]);
			print("distEval", 'INDEL', '%FNR', (100 * FN[1] / (FN[1] + TP[1])).toFixed(2));
			print("distEval", 'INDEL', '%FDR', (100 * FP[1] / (FP[1] + TP1[1])).toFixed(2));
		}
	}
	if (hp != null) hp.destroy();
}

/***********************
 *** Main() function ***
 ***********************/

function main(args)
{
	if (args.length == 0) {
		print("\nUsage:    k8 hapdip.js <command> [arguments]");
		print("Version:  r56\n");
		print("Commands: eval     evaluate a pair of CHM1 and NA12878 VCFs");
		print("          distEval distance-based VCF comparison");
		print("");
		print("          deovlp   remove overlaps between variants");
		print("          upd1gt   update genotypes in a single-sample VCF");
		print("          anno     recompute some INFO");
		print("          filter   filter anno'd VCF");
		print("");
		print("          qst1     vcf stats stratified by QUAL, one sample only");
		print("          vcf2bed  convert VCF to unary BED");
		print("          ucf2mat  convert unary VCF to bit matrix");
		print("          vcfsub   subset, reorder and rename samples in VCF");
		print("          vcfsum   write total allele/genotype counts and depths in INFO");
		print("          vcfswap  swap the two alleles in an unphased genotype");
		print("");
		print("          cg2vcf   convert CG's masterVarBeta to VCF");
		print("          bedovlp  count lines overlapping in a second bed");
		print("          bedcmpm  compare multiple sorted BED files");
		print("          cbs      circular binary segmentation");
		print("");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'eval') b8_eval(args);
	else if (cmd == 'deovlp') b8_deovlp(args);
	else if (cmd == 'upd1gt') b8_upd1gt(args);
	else if (cmd == 'anno') b8_anno(args);
	else if (cmd == 'filter') b8_filter(args);
	else if (cmd == 'qst1') b8_qst1(args);
	else if (cmd == 'cg2vcf') b8_cg2vcf(args);
	else if (cmd == 'bedovlp') b8_bedovlp(args);
	else if (cmd == 'bedcmpm') b8_bedcmpm(args);
	else if (cmd == 'cbs') b8_cbs(args);
	else if (cmd == 'distEval') b8_distEval(args);
	else if (cmd == 'vcfsub') b8_vcfsub(args);
	else if (cmd == 'vcfsum') b8_vcfsum(args);
	else if (cmd == 'vcf2bed') b8_vcf2bed(args);
	else if (cmd == 'vcfswap') b8_vcfswap(args);
	else if (cmd == 'postatom' || cmd == 'atompost') b8_postatom(args);
	else if (cmd == 'atomcnt') b8_atomcnt(args);
	else warn("Unrecognized command");
}

main(arguments);
