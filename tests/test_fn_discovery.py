"""FN discovery tests: targeted locus review, filtered recovery, classification."""
from avrs.fn_discovery import (load_hotspot_file, targeted_locus_review,
                               recover_filtered_candidates, summarize)
from avrs.db.models import FNDiscoveryCandidate
from avrs.io.fasta import Fasta
from avrs.variants import import_vcf
from avrs.review import import_sample_bam_vcf


def test_hotspot_loader_reads_tsv(tmp_path):
    p = tmp_path / "hs.tsv"
    p.write_text("#chrom\tpos\tref\talt\nchrDemo1\t20000\tG\tT\n", encoding="utf-8")
    loci = load_hotspot_file(p)
    assert len(loci) == 1
    assert loci[0].chrom == "chrDemo1" and loci[0].pos == 20000
    assert loci[0].ref == "G" and loci[0].alt == "T"


def test_targeted_locus_finds_caller_miss(cfg, db_session, demo_dataset):
    """Scenario 6 is present in the BAM but absent from the VCF."""
    sample_id = demo_dataset["sample_id"]
    smp = import_sample_bam_vcf(db_session, cfg, sample_id, demo_dataset["bam"],
                                demo_dataset["vcf"])
    fasta = Fasta(demo_dataset["fasta"])
    import_vcf(db_session, demo_dataset["vcf"], smp.id, sample_id,
               demo_dataset["reference_build"], fasta=fasta)
    db_session.commit()

    loci = load_hotspot_file(demo_dataset["hotspots"])
    cands = targeted_locus_review(db_session, sample_id, cfg, demo_dataset["bam"],
                                  demo_dataset["fasta"], loci)
    assert cands
    miss = [c for c in cands if c.possible_missed]
    assert miss, "the caller-miss scenario must surface as a candidate"
    assert miss[0].classification == "possible_missed_variant"
    # must NOT be reported as a confirmed false negative
    s = summarize(db_session, sample_id)
    assert s["possible_missed"] >= 1
    assert "NOT a confirmed false negative" in s["disclaimer"]
    rows = db_session.query(FNDiscoveryCandidate).all()
    assert any(r.possible_missed_variant for r in rows)


def test_filtered_recovery_creates_candidates(cfg, db_session, demo_dataset):
    sample_id = demo_dataset["sample_id"]
    smp = import_sample_bam_vcf(db_session, cfg, sample_id, demo_dataset["bam"],
                                demo_dataset["vcf"])
    cands = recover_filtered_candidates(db_session, sample_id, cfg,
                                        demo_dataset["vcf"], demo_dataset["bam"],
                                        demo_dataset["fasta"])
    assert cands
    # a FILTER-failed call is a candidate, never auto-restored to PASS
    for c in cands:
        assert c.source == "filtered_recovery"
        assert c.caller_status == "filtered_in_raw"
    rows = db_session.query(FNDiscoveryCandidate).all()
    assert all(r.possible_missed_variant in (True, False) for r in rows)


def test_in_classified_locus_is_callable_and_in_vcf(cfg, db_session, demo_dataset):
    sample_id = demo_dataset["sample_id"]
    smp = import_sample_bam_vcf(
        db_session, cfg, sample_id, demo_dataset["bam"], demo_dataset["vcf"])
    fasta = Fasta(demo_dataset["fasta"])
    import_vcf(db_session, demo_dataset["vcf"], smp.id, sample_id,
               demo_dataset["reference_build"], fasta=fasta)
    db_session.commit()
    from avrs.fn_discovery import LocusSpec
    from avrs.db.models import Variant
    v = (db_session.query(Variant)
         .filter_by(sample_id=smp.id, chrom="chrDemo1", pos=5000).one())
    loci = [LocusSpec(chrom="chrDemo1", pos=5000, ref=v.ref, alt=v.alt)]
    cands = targeted_locus_review(db_session, sample_id, cfg, demo_dataset["bam"],
                                  demo_dataset["fasta"], loci)
    assert cands[0].classification == "in_vcf"
    assert cands[0].possible_missed is False
