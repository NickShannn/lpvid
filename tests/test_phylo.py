from cdc_pipeline.stages.phylo import phylo_cache_key


def test_cache_key_changes_when_adapter_content_changes(tmp_path):
    base = tmp_path / "base"
    adapter = tmp_path / "adapter"
    base.mkdir()
    adapter.mkdir()
    (base / "model.safetensors").write_bytes(b"base")
    (base / "config.json").write_text("{}")
    (adapter / "adapter_model.safetensors").write_bytes(b"v1")
    first = phylo_cache_key(base, adapter)
    (adapter / "adapter_model.safetensors").write_bytes(b"v2")
    second = phylo_cache_key(base, adapter)
    assert first != second
