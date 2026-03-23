"""Tests for the bam_to_bed AWK-based SAM-to-BED conversion."""

import pytest

from oligominer.specificity.alignment.bam_to_bed import sam_to_bed


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _make_sam_line(name, flag, chrom, pos, cigar, score=0):
    """Build a minimal SAM record with the fields the AWK script reads."""
    seq = "A" * 20
    qual = "*"
    return "\t".join([
        name, str(flag), chrom, str(pos), "42", cigar,
        "*", "0", "0", seq, qual, "AS:i:%d" % score,
    ])


def _parse_bed_line(line):
    """Return a dict of the BED fields produced by the AWK script."""
    parts = line.split("\t")
    return {
        "chrom": parts[0],
        "start": int(parts[1]),
        "end": int(parts[2]),
        "name": parts[3],
        "score": parts[4],
        "strand": parts[5],
        "cigar": parts[6],
    }


# ---------------------------------------------------------------------------
# strand extraction (the bug-fix under test)
# ---------------------------------------------------------------------------

class TestStrandExtraction:
    """Verify the portable bit-check for the 0x10 reverse-complement flag."""

    @pytest.mark.parametrize("flag, expected_strand", [
        (0, "+"),       # no flags
        (1, "+"),       # paired
        (2, "+"),       # proper pair
        (4, "+"),       # unmapped
        (16, "-"),      # reverse complement
        (17, "-"),      # paired + reverse
        (32, "+"),      # mate reverse (not self)
        (48, "-"),      # mate reverse + self reverse
        (99, "+"),      # typical forward read in proper pair
        (147, "-"),     # typical reverse read in proper pair
        (163, "+"),     # another common forward flag
        (256, "+"),     # secondary alignment, forward
        (272, "-"),     # secondary alignment, reverse
    ])
    def test_strand_from_flag(self, flag, expected_strand):
        sam = _make_sam_line("probe", flag, "chr1", 100, "20M") + "\n"
        result = sam_to_bed(sam)
        bed = _parse_bed_line(result.strip())
        assert bed["strand"] == expected_strand


# ---------------------------------------------------------------------------
# coordinate handling
# ---------------------------------------------------------------------------

class TestCoordinates:

    def test_simple_match(self):
        """20M cigar at pos 100 → BED [99, 119)."""
        sam = _make_sam_line("p1", 0, "chr1", 100, "20M") + "\n"
        bed = _parse_bed_line(sam_to_bed(sam).strip())
        assert bed["start"] == 99
        assert bed["end"] == 119

    def test_soft_clip_adjusts_start(self):
        """5S15M at pos 200 → start backed up by 5."""
        sam = _make_sam_line("p2", 0, "chr1", 200, "5S15M") + "\n"
        bed = _parse_bed_line(sam_to_bed(sam).strip())
        assert bed["start"] == 194   # (200 - 5) - 1 = 194

    def test_start_clamp_at_zero(self):
        """Soft-clip that would push start below 0 is clamped."""
        sam = _make_sam_line("p3", 0, "chr1", 3, "10S10M") + "\n"
        bed = _parse_bed_line(sam_to_bed(sam).strip())
        assert bed["start"] == 0


# ---------------------------------------------------------------------------
# score parsing
# ---------------------------------------------------------------------------

class TestScoreParsing:

    def test_alignment_score(self):
        sam = _make_sam_line("p1", 0, "chr1", 100, "20M", score=-5) + "\n"
        bed = _parse_bed_line(sam_to_bed(sam).strip())
        assert bed["score"] == "-5"

    def test_zero_score(self):
        sam = _make_sam_line("p1", 0, "chr1", 100, "20M", score=0) + "\n"
        bed = _parse_bed_line(sam_to_bed(sam).strip())
        assert bed["score"] == "0"


# ---------------------------------------------------------------------------
# multi-record input
# ---------------------------------------------------------------------------

class TestMultiRecord:

    def test_two_records(self):
        lines = "\n".join([
            _make_sam_line("fwd", 0, "chr1", 100, "20M"),
            _make_sam_line("rev", 16, "chr2", 200, "5S15M", score=-3),
        ]) + "\n"
        beds = [_parse_bed_line(l) for l in sam_to_bed(lines).strip().split("\n")]
        assert len(beds) == 2
        assert beds[0]["strand"] == "+"
        assert beds[1]["strand"] == "-"
        assert beds[0]["chrom"] == "chr1"
        assert beds[1]["chrom"] == "chr2"
