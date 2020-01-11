from quantized.config import to_bool, show, conf, Config
from quantized import config

from pytest import raises


def test_to_bool():

    trues = [True, "True", "true", "yes", "Yes", "YES", "TRUE"]
    falses = [False, "False", "notrue", "No"]

    for v in trues:
        assert to_bool(v)

    for v in falses:
        assert not to_bool(v)

    with raises(TypeError):
        to_bool(["list"])


def test_show():
    show(conf)
    show()


def test_set_conf(tmp_path):
    p = tmp_path / "conf.yaml"

    c = Config(harmonic_oscillator_max_n=20_000_000)
    config.set_conf(c, path=p)
    loaded = config.load(p)
    assert loaded == c


def test_restore_defaults(tmp_path):

    p = tmp_path / "conf.yaml"

    n = 2_334_121_123
    config.set_items(p, harmonic_oscillator_max_n=n)
    c = config.load(p)
    assert c.harmonic_oscillator_max_n == n

    config.restore_defaults(p)
    reloaded = config.load(p)
    assert reloaded == Config()


def test_set_items(tmp_path):
    p = tmp_path / "conf.yaml"

    n = 25_000_000
    config.set_items(p, harmonic_oscillator_max_n=n)
    new_conf = config.load(p)
    assert new_conf.harmonic_oscillator_max_n == n


def test_clear_cache(tmp_path):
    cache_path = tmp_path / "cache"
    cache_path.mkdir()

    assert cache_path.is_dir()

    c = Config(cache_dir=cache_path)
    config.clear_cache(c)
    assert not cache_path.is_dir()
