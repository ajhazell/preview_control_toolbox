function tests = testPCTDemo
    tests = functiontests(localfunctions);
end

% TODO break the checks in PCTDemo into individual tests

function tests = testBasicCheckNoErrors(testCase)
% basic check that no errors are thrown

    run('PCTDemo')
end