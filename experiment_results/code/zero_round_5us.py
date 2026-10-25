"""Same logical X/Z zero-round experiment, with 5 us storage.

python code/zero_round_5us.py --self-test
python code/zero_round_5us.py --prepare
python code/zero_round_5us.py --hardware --resume-dir DIR
python code/reproduce_analysis.py

Reuses preparation, final readout, mapping, calibration and decoder utilities
from the near-zero experiment. Only --hardware submits real circuits.
"""

from pathlib import Path
import sys

from zero_round import main


if __name__ == '__main__':
    main(storage_ns=5000., results_folder='zero_round_5us')
