"""Test all / specified contributed problems, and this depends on:
1. results generated by: run_examples.py
2. criteria specified in: criteria.py
"""

import pathlib
import pytest

import run_examples
import validate
import report


args = validate.args

def _pre_build():
    return args.pre or (not args.pre and not args.post)

@pytest.fixture
def pre_build():
    return _pre_build()

def _all_contribs():
    pre = _pre_build()
    problems = run_examples.problems_to_run(args.contribs)
    print("🥃 Running " + ("pre-" if pre else "post-") + "build tests for the following contributions:")
    print("- " + "\n- ".join([c[0] for c in problems]) + "\n")
    # results = run_examples.run_problems(problems, pre_build=pre)
    # return results
    return problems

@pytest.fixture(params=_all_contribs())
def contrib(request):
    return request.param

def test_contrib(contrib, pre_build):
    _report = report.compliance_report([contrib[0]], pre_build)
    for _r in _report.values():
        report.pprint_compliance_report(_report)
        if isinstance(_r["api_compliance"], Exception):
            raise _r["api_compliance"]
        assert _r["api_compliance"], \
            "Not API-compliant. Check report above for details."

def main():
    pytest.main([pathlib.Path(__file__)])

if __name__ == "__main__":
    main()
