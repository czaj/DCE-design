% EXAMPLE_1_BASIC  Efficient design for a multinomial logit model with fixed priors.
%
% Two unlabelled hotels described by price, stars, meals and pool, plus the
% purpose of the trip as a scenario variable that changes the price effect.
% 18 choice tasks in 3 blocks. Takes a few seconds.

addpath(fileparts(fileparts(mfilename('fullpath'))));

spec = struct();
spec.name = 'hotel';
spec.alternatives = {'hotel1', 'hotel2'};
spec.rows = 18;
spec.blocks = 3;

spec.attributes.purpose = [0 1];             % 0 leisure, 1 business
spec.attributes.price   = [240 180 120];
spec.attributes.stars   = [1 3 5];           % base level of a coded attribute: the first listed
spec.attributes.meals   = [0 1 2];
spec.attributes.pool    = [0 1];
spec.scenario = {'purpose'};                 % one column, the same for both hotels
spec.labels.meals = {'none', 'breakfast', 'all-inclusive'};

u = 'b_price*price + b_stars.dummy*stars + b_meals.dummy*meals + b_pool.dummy*pool + b_pp*price*purpose';
spec.model.utility.hotel1 = u;
spec.model.utility.hotel2 = u;
spec.model.priors.b_price = -0.02;
spec.model.priors.b_stars = [0.4 0.7];       % 3 and 5 stars against 1 star
spec.model.priors.b_meals = [0.3 0.8];       % breakfast and all-inclusive against none
spec.model.priors.b_pool  = 0.3;
% b_pp has no prior: it is set to zero and the result says so

dce_check(spec);                              % how the specification is read
result = dce_design(spec);                  % search, block, evaluate, print the summary

disp(result.design);
% dce_export(result, fullfile(tempdir, 'hotel_design'));   % csv, xlsx, json, mat, report
