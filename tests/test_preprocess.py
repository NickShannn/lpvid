from cdc_pipeline.stages.preprocess import fastp_args, megahit_args, rrna_args


CFG = {
    "qc": {"threads": 16},
    "assembly": {"threads": 64, "min_contig_len": 600, "memory_gb": 64},
}


def test_megahit_supports_se_and_pe_without_fabricating_r2(tmp_path):
    se = megahit_args(["r1.fq"], CFG, tmp_path / "se")
    pe = megahit_args(["r1.fq", "r2.fq"], CFG, tmp_path / "pe")
    assert se[-2:] == ["-r", "r1.fq"]
    assert pe[-4:] == ["-1", "r1.fq", "-2", "r2.fq"]
    assert "--min-contig-len" in se
    assert "600" in se


def test_fastp_and_rrna_commands_have_explicit_se_pe_contract(tmp_path):
    se = fastp_args("r1.fq", None, CFG, tmp_path)
    pe = fastp_args("r1.fq", "r2.fq", CFG, tmp_path)
    assert "-I" not in se
    assert pe[pe.index("-I") + 1] == "r2.fq"
    rrna_se = rrna_args("index", "r1.fq", None, tmp_path, "ssu")
    rrna_pe = rrna_args("index", "r1.fq", "r2.fq", tmp_path, "ssu")
    assert rrna_se[-2:] == ["-U", "r1.fq"]
    assert "--un-gz" in rrna_se
    assert rrna_pe[-4:] == ["-1", "r1.fq", "-2", "r2.fq"]
    assert "--un-conc-gz" in rrna_pe
