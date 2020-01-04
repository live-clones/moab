import sys
import traceback

if sys.version_info < (3, 0):
    from collections import Iterable
else:
    from collections.abc import Iterable

class colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def test_driver(test_list):
    ret_val = 0
    for test in test_list:
        try:
            test()
        except:
            print(colors.FAIL + "FAIL" + colors.ENDC + ": " + test.__name__)
            ret_val += 1
            traceback.print_exc()
        else:
            print(colors.OKGREEN + "PASS" + colors.ENDC + ": " + test.__name__)
    sys.exit(ret_val)

def CHECK_ITER_EQ(actual_value, expected_value):
    CHECK_EQ(len(actual_value), len(expected_value))
    for a,e in zip(actual_value, expected_value):
        if isinstance(a, str) and (e, str):
            CHECK_EQ(a,e)
            continue
        if isinstance(a, Iterable) and isinstance(e, Iterable):
            CHECK_ITER_EQ(a,e)
        else:
            CHECK_EQ(a,e)

def CHECK_EQ(actual_value, expected_value):
    err_msg = "Expected value: {} Actual value: {}"
    err_msg = err_msg.format(expected_value, actual_value)
    result = expected_value == actual_value
    assert result, err_msg

def CHECK(actual_value):
    CHECK_EQ(actual_value, True)

def CHECK_NOT(actual_value):
    CHECK_EQ(actual_value, False)

def CHECK_NOT_EQ(actual_value, expected_value):
    err_msg = "Expected value: not {} Actual value: {}"
    err_msg = err_msg.format(expected_value, actual_value)
    result = expected_value != actual_value
    assert result, err_msg

def CHECK_TYPE(actual_type, expected_type):
    err_msg = "Type {} does not match expected type {}"
    err_msg = err_msg.format(type(actual_type), expected_type)
    result = type(actual_type) == expected_type
    assert result, err_msg
