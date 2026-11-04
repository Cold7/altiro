import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from orthoresponse import classify_responses  # noqa: E402
from orthoresponse.plotting import plot_model_patterns  # noqa: E402


def _ticks(fig):
    ax = next(a for a in fig.axes if a.get_xticklabels() and a.get_xticklabels()[0].get_text())
    return [t.get_text() for t in ax.get_xticklabels()]


def test_default_condition_labels():
    fig = plot_model_patterns()
    assert _ticks(fig) == ["C", "S", "C", "S"]
    assert fig._suptitle.get_text() == "Expression models (C = control, S = stress)"


def test_custom_condition_labels_and_title(contrasts):
    cls = classify_responses(contrasts, group_names=("Sly", "Spenn"))
    fig = plot_model_patterns(cls, condition_labels=("Ctrl", "Drought"))
    assert _ticks(fig) == ["Ctrl", "Drought", "Ctrl", "Drought"]
    assert fig._suptitle.get_text().startswith("Expression models (Ctrl = control, Drought = stress)")
    assert [t.get_text() for t in fig.legends[0].get_texts()] == ["Sly", "Spenn"]
    fig = plot_model_patterns(cls, condition_labels=["C", "D"], title="Tomato drought")
    assert fig._suptitle.get_text() == "Tomato drought"


def test_condition_labels_validation():
    with pytest.raises(ValueError):
        plot_model_patterns(condition_labels=("Control",))
