import gzip
import logging
from collections import defaultdict
from typing import Optional, DefaultDict

import pyfaidx
from dataclasses import dataclass


class FastaNotIndexedError(Exception):
    pass


class InvalidRegion(Exception):
    pass

class NumericSampleIds:
	"""
	Mapping of sample names (strings) to numeric ids.
	"""
	def __init__(self):
		self.mapping = {}
		self.frozen = False

	def __getitem__(self, sample):
		if not self.frozen and sample not in self.mapping:
			self.mapping[sample] = len(self.mapping)
		return self.mapping[sample]

	def __len__(self):
		return len(self.mapping)

	def __str__(self):
		return str(self.mapping)

	def freeze(self):
		"""No longer allow modifications"""
		# TODO We should try to have a second class (FrozenNumericSampleIds) for
		# this or try to get rid of NumericSampleIds
		self.frozen = True

	def inverse_mapping(self):
		"""Returns a dict mapping numeric ids to sample names."""
		return {numeric_id:name for name, numeric_id in self.mapping.items()}
	
	def __getstate__(self):
		return (self.mapping, self.frozen)
	
	def __setstate__(self, state):
		mapping, frozen = state
		self.mapping = mapping
		self.frozen = frozen

def detect_file_format(path):
    """
    Detect file format and return 'BAM', 'CRAM', 'VCF' or None. None indicates an
    unrecognized file format.

    'VCF' is returned for both uncompressed and compressed VCFs (.vcf and .vcf.gz).
    """
    with open(path, "rb") as f:
        first_bytes = f.read(16)
        if first_bytes.startswith(b"CRAM"):
            return "CRAM"
        if first_bytes.startswith(b"##fileformat=VCF"):
            return "VCF"

    gzip_header = b"\037\213"
    if first_bytes.startswith(gzip_header):
        with gzip.GzipFile(path, "rb") as f:
            first_bytes = f.read(16)
            if first_bytes.startswith(b"BAM\1"):
                return "BAM"
            elif first_bytes.startswith(b"##fileformat=VCF"):
                return "VCF"

    return None


def IndexedFasta(path):
    try:
        f = pyfaidx.Fasta(path, as_raw=True, sequence_always_upper=True, build_index=False)
    except pyfaidx.IndexNotFoundError:
        raise FastaNotIndexedError(path)
    return f


def plural_s(n: int) -> str:
    return "" if n == 1 else "s"


@dataclass
class Region:
    chromosome: str
    start: int
    end: Optional[int]

    def __repr__(self):
        return f'Region("{self.chromosome}", {self.start}, {self.end})'

    @staticmethod
    def parse(spec: str):
        """
        >>> Region.parse("chr1")
        Region("chr1", 0, None)
        >>> Region.parse("chr1:")
        Region("chr1", 0, None)
        >>> Region.parse("chr1:101")
        Region("chr1", 100, None)
        >>> Region.parse("chr1:101-")
        Region("chr1", 100, None)
        >>> Region.parse("chr1:101-200")
        Region("chr1", 100, 200)
        >>> Region.parse("chr1:101:200")  # for backwards compatibility
        Region("chr1", 100, 200)
        """
        parts = spec.split(":", maxsplit=1)
        chromosome = parts[0]
        if len(parts) == 1 or not parts[1]:
            start, end = 0, None
        else:
            try:
                sep = ":" if ":" in parts[1] else "-"
                start_end = parts[1].split(sep, maxsplit=1)
                start = int(start_end[0]) - 1
                if len(start_end) == 1 or not start_end[1]:
                    end = None
                else:
                    end = int(start_end[1])
                    if end <= start:
                        raise InvalidRegion("end is before start in specified region")
            except ValueError:
                raise InvalidRegion("Region must be specified as chrom[:start[-end]])") from None
        return Region(chromosome, start, end)


_warning_count: DefaultDict[str, int] = defaultdict(int)


def warn_once(logger, msg: str, *args) -> None:
    if _warning_count[msg] == 0 and not logger.isEnabledFor(logging.DEBUG):
        logger.warning(msg + " Hiding further warnings of this type, use --debug to show", *args)
    else:
        logger.debug(msg, *args)
    _warning_count[msg] += 1
