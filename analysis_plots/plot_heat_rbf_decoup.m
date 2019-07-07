clc; close all; clear all;

%% problem properties
% % problem properties
% pp = containers.Map({            'prblm', 'dof per node', 'conductivity'}, ...
%                      { "steady heat MMS",              1,              1});
%                  
% % background mesh
% box1 = [-2,2; -2,2];
% h1   = [ 0.2,  0.2]/2/2/2;
% 
% % boundary condition map - bottom, right, top, left
% bc1 = containers.Map({   'bottom',     'right',       'top',      'left'}, ...
%                      {"dirichlet", "dirichlet", "dirichlet", "dirichlet"});
% % near body mesh
% box2 = [-1.13625,0.86375; -1.13625,0.86375];
% h2   = [             0.2,              0.2]/2/2/2;
% 
% % boundary condition map - bottom, right, top, left
% bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
%                      {"overset", "overset", "overset", "overset"});
%
% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order', 'solve type'}, ...
%                          { 2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 2, "decoupled" });
% 
% % convergence parameters
% N_iters    =    10; % Maximum number of Newton steps
% resnrmdrop = 1e-09; % Newton convergence criteria

del    = 0.2;
h      = [del del/2 del/4 del/8];
L2_p2  = [5.152661e-02 9.918104e-03 2.935698e-03 7.070890e-04];
L2_p2_1em4  = [0.051548827215790 0.009916565074277 0.002945711425439 8.744870950362990e-04];
L2_p2_1em5  = [5.1529189759193503e-02, 9.9178520555716614e-03, 2.9366559655353151e-03, 7.2237730065284026e-04];
y_vals2= 0.1*h.^2;

figure(1)
clf
hold on
plot(h, L2_p2,'--*r','LineWidth', 2);
plot(h, L2_p2_1em4,'--*b','LineWidth', 2);
plot(h, L2_p2_1em5,'--*c','LineWidth', 2);
plot(h,y_vals2,'-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('RBF$(C^2)$, coupled','RBF$(C^2)$, tol=1e-04','RBF$(C^2)$, tol=1e-05','slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

L2_C2_hp2 = [0.111185401437300,0.0686941869452710,0.0681288190682004,0.0560122111260676,0.0566958407943873,0.0528484946508973,0.0532169058847265,0.0519510325673926,0.0520934570342101,0.0516680991851441,0.0517186119779739,0.0515744421926258,0.0515918796179365,0.0515428613235235,0.0515488272157900,0.0515321425664045,0.0515341775302890,0.0515284963450578,0.0515291897591935,0.0515272550465859,0.0515274912451104,0.0515268323550703,0.0515269128021149,0.0515266884057147,0.0515267158041300,0.0515266393816671,0.0515266487128094,0.0515266226856410,0.0515266258635550,0.0515266169994872,0.0515266180817900,0.0515266150629554,0.0515266154315553,0.0515266144034311,0.0515266145289652,0.0515266141788170,0.0515266142215701,0.0515266141023202,0.0515266141168806,0.0515266140762677,0.0515266140812265,0.0515266140673950,0.0515266140690838,0.0515266140643732,0.0515266140649483,0.0515266140633441,0.0515266140635400,0.0515266140629936,0.0515266140630603,0.0515266140628742,0.0515266140628969,0.0515266140628336,0.0515266140628413,0.0515266140628197,0.0515266140628224,0.0515266140628150,0.0515266140628159,0.0515266140628134,0.0515266140628137,0.0515266140628128,0.0515266140628129,0.0515266140628127,0.0515266140628127,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628125,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628125,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628125,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126,0.0515266140628126];
L2_C2_hp1 = [0.0930988430622235,0.0665501374817493,0.0531333037032679,0.0429427610956081,0.0345882614555620,0.0286120419337310,0.0233720839949760,0.0199071573988117,0.0167258159003750,0.0148570779559416,0.0130527004846358,0.0121610129736808,0.0112204790300622,0.0108585555050914,0.0104000538864896,0.0102829984156859,0.0100661234491174,0.0100456974419397,0.00994280186585982,0.00995390024553611,0.00990383228468514,0.00992139114297127,0.00989607622745086,0.00991182560157430,0.00989844503080072,0.00991048998845947,0.00990309470750287,0.00991167490266841,0.00990741918022674,0.00991330675509680,0.00991077420747525,0.00991472846016742,0.00991318166933887,0.00991580346640394,0.00991484051902438,0.00991656507427656,0.00991595744779656,0.00991708615473643,0.00991669916423271,0.00991743554565469,0.00991718752923258,0.00991766697458819,0.00991750736159707,0.00991781911116955,0.00991771610849410,0.00991791864586076,0.00991785205557166,0.00991798356803968,0.00991794046733796,0.00991802583138082,0.00991799791295845,0.00991805330960270,0.00991803521645381,0.00991807116053764,0.00991805943106910,0.00991808275112318,0.00991807514551345,0.00991809027431651,0.00991808534202018,0.00991809515638050,0.00991809195746434,0.00991809832407106,0.00991809624924552,0.00991810037921317,0.00991809903342545,0.00991810171247335,0.00991810083953803,0.00991810257738366,0.00991810201115167,0.00991810313845267,0.00991810277116089,0.00991810350241327,0.00991810326416438,0.00991810373850884,0.00991810358396476,0.00991810389165931,0.00991810379141146,0.00991810399100451,0.00991810392597679,0.00991810405544728,0.00991810401326571,0.00991810409724963,0.00991810406988769,0.00991810412436568,0.00991810410661679,0.00991810414195514,0.00991810413044193,0.00991810415336491,0.00991810414589663,0.00991810416076613,0.00991810415592164,0.00991810416556709,0.00991810416242462,0.00991810416868133,0.00991810416664291,0.00991810417070146,0.00991810416937919,0.00991810417201186,0.00991810417115415,0.00991810417286189,0.00991810417230550,0.00991810417341326,0.00991810417305235,0.00991810417377093,0.00991810417353681,0.00991810417400293,0.00991810417385107,0.00991810417415344,0.00991810417405493,0.00991810417425105,0.00991810417418717,0.00991810417431438,0.00991810417427293,0.00991810417435547,0.00991810417432858,0.00991810417438211,0.00991810417436467,0.00991810417439940,0.00991810417438808,0.00991810417441061,0.00991810417440327,0.00991810417441787,0.00991810417441313,0.00991810417442259,0.00991810417441950,0.00991810417442565,0.00991810417442365,0.00991810417442763,0.00991810417442635,0.00991810417442893,0.00991810417442808,0.00991810417442976,0.00991810417442921,0.00991810417443031,0.00991810417442994,0.00991810417443067,0.00991810417443042,0.00991810417443089,0.00991810417443073,0.00991810417443104,0.00991810417443093,0.00991810417443113,0.00991810417443107,0.00991810417443118,0.00991810417443115,0.00991810417443124,0.00991810417443120,0.00991810417443125,0.00991810417443123,0.00991810417443127,0.00991810417443126,0.00991810417443128,0.00991810417443127,0.00991810417443129,0.00991810417443129,0.00991810417443129,0.00991810417443130,0.00991810417443129,0.00991810417443130,0.00991810417443129];
L2_C2_hp05 = [0.108872378292197,0.0716581409843094,0.0731566866700331,0.0509826853166998,0.0542561028522759,0.0381534809230521,0.0415851252455452,0.0293066534787208,0.0325644159896145,0.0229811093455545,0.0259091936825891,0.0183239049220833,0.0208723245403704,0.0148093337550688,0.0169897952272323,0.0121075975482713,0.0139593948354851,0.0100059825812299,0.0115756386093749,0.00836178641461778,0.00969311000927417,0.00707479755304768,0.00820513313485946,0.00607126416064889,0.00703095600060575,0.00529444185952012,0.00610777680105078,0.00469895550953250,0.00538563026550981,0.00424748519232337,0.00482405585021889,0.00390892990102055,0.00438993250477700,0.00365746239348537,0.00405607983617497,0.00347200483764174,0.00380032465634694,0.00333578853724760,0.00360480401565374,0.00323583589267339,0.00345535773618969,0.00316235344864808,0.00334094619583215,0.00310810018412533,0.00325308890524217,0.00306780123049459,0.00318534585660252,0.00303765014641285,0.00313286419848405,0.00301491299139857,0.00309200174346384,0.00299762756296105,0.00306002700037976,0.00298438208121382,0.00303488756704440,0.00297415576220755,0.00301503526905246,0.00296620550824884,0.00299929609740863,0.00295998598547464,0.00298677432195707,0.00295509342899212,0.00297678208623293,0.00295122612001898,0.00296878772621020,0.00294815650501816,0.00296237773328088,0.00294571142543920,0.00295722862452109,0.00294375800168161,0.00295308600766645,0.00294219347104893,0.00294974888641079,0.00294093780479035,0.00294705780210295,0.00293992829158875,0.00294488580285390,0.00293911552369938,0.00294313151259074,0.00293846039268661,0.00294171377278956,0.00293793181898350,0.00294056747208575,0.00293750502028664,0.00293964028067997,0.00293716017970491,0.00293889007947132,0.00293688141349848,0.00293828292661989,0.00293665596553532,0.00293779144269260,0.00293647357489002,0.00293739352381377,0.00293632597678339,0.00293707131321800,0.00293620650699017,0.00293681037731732,0.00293610978707573,0.00293659904427502,0.00293603147314844,0.00293642787214332,0.00293596805477612,0.00293628922059178,0.00293591669369482,0.00293617690565767,0.00293587509419787,0.00293608592116367,0.00293584139882565,0.00293601221375775,0.00293581410431293,0.00293595250114046,0.00293579199379062,0.00293590412511368,0.00293577408205346,0.00293586493273114,0.00293575957134520,0.00293583318014402,0.00293574781562079,0.00293580745478680,0.00293573829164699,0.00293578661239100,0.00293573057562461,0.00293576972599203,0.00293572432427164,0.00293575604463958,0.00293571925951216,0.00293574495996096,0.00293571515608165,0.00293573597908139,0.00293571183149133,0.00293572870269067,0.00293570913790154,0.00293572280727732,0.00293570695554044,0.00293571803073740,0.00293570518737354,0.00293571416071668,0.00293570375478621,0.00293571102516646,0.00293570259408636,0.00293570848469273,0.00293570165367137,0.00293570642635779,0.00293570089173307,0.00293570475865869,0.00293570027439852,0.00293570340745880,0.00293569977422376,0.00293570231269176,0.00293569936897349,0.00293570142569100,0.00293569904063249,0.00293570070702614,0.00293569877460460,0.00293570012475009,0.00293569855906385,0.00293569965297876,0.00293569838442867,0.00293569927074043,0.00293569824293591];
L2_C2_hp025 = [0.116160955553106,0.0966550919124544,0.0962366682772385,0.0817324949838165,0.0838119009650016,0.0711097663591314,0.0737838821608827,0.0623420607678322,0.0652396865476925,0.0548766015708811,0.0578541622468741,0.0484596836964612,0.0514248412216742,0.0429124392298928,0.0457990460560758,0.0380932629273440,0.0408556019851734,0.0338875924478768,0.0364964170252739,0.0302022612746206,0.0326410473274240,0.0269611533117610,0.0292227930265073,0.0241016992958935,0.0261857994814083,0.0215720887683794,0.0234828709473100,0.0193290822804060,0.0210738020819901,0.0173363069286102,0.0189240910517496,0.0155629306159691,0.0170039368756133,0.0139826287571923,0.0152874504828244,0.0125727753222139,0.0137520278286224,0.0113138056682732,0.0123778468806945,0.0101887110703607,0.0111474600250764,0.00918263451952756,0.0100454605540336,0.00828254472736100,0.00905820713893877,0.00747697085133380,0.00817359407604410,0.00675578466004332,0.00738085798719707,0.00611002002270057,0.00667041382249085,0.00553172199019133,0.00603371463736410,0.00501381953050998,0.00546313084177762,0.00455001733637362,0.00495184554623266,0.00413470314787552,0.00449376333439235,0.00376286781107004,0.00408343033124121,0.00343003588656674,0.00371596385099634,0.00313220507710708,0.00338699023165102,0.00286579309451739,0.00309258971643962,0.00262759086029124,0.00282924744407552,0.00241472114951139,0.00259380977222760,0.00222460195881838,0.00238344529175822,0.00205491401557563,0.00219560999936307,0.00190357195405625,0.00202801618778533,0.00176869876982225,0.00187860468826525,0.00164860322833500,0.00174552016048839,0.00154175995022135,0.00162708917129747,0.00144679192541677,0.00152180083486851,0.00136245522404427,0.00142828980438456,0.00128762567674266,0.00134532140998919,0.00122128729560586,0.00127177873294457,0.00116252220391060,0.00120665139585493,0.00111050184322319,0.00114902583885288,0.00106447923405880,0.00109807684703165,0.00102378208297278,0.00105306009921219,0.000987806554588361,0.00101330552432069,0.000956011559430462,0.000978211278673489,0.000927913444067539,0.000947238192514025,0.000903081005093034,0.000919904572986622,0.000881130779515555,0.000895781288781470,0.000861722588855944,0.000874487095036299,0.000844555331657343,0.000855684183155743,0.000829363029383758,0.000839073957984927,0.000815911134863406,0.000824393054598749,0.000803993112039508,0.000811409610150002,0.000793427292435608,0.000799919804543030,0.000784054008872398,0.000789744679004583,0.000775733001724978,0.000780727235518929,0.000768341088168337,0.000772729813809845,0.000761770080881977,0.000765631736906552,0.000755924939768601,0.000759327211776838,0.000750722138396461,0.000753723468216107,0.000746088225994288,0.000748739117115912,0.000741958565745858,0.000744302708254163,0.000738276230664336,0.000740351467649107,0.000734991039297253,0.000736830195086808,0.000732058714767017,0.000733690303468100,0.000729440152062061,0.000730888982954290,0.000727100779962478,0.000728388474383993,0.000725010005439724,0.000726155437981880,0.000723140729763697,0.000724160404907940,0.000721468926850698,0.000722377300652840,0.000719973275576286,0.000720783030639200,0.000718634838849310,0.000719357119623457,0.000717436783199028,0.000718081397602181,0.000716364133469764];

iterations = 1:length(L2_C2_hp025);
oneD_arr = ones(length(L2_C2_hp025));
y_vals = h.^2;

figure(2)
clf
h1 = plot3(iterations,0.2*oneD_arr,L2_C2_hp2,'--*r','LineWidth', 2);
hold on
h2 = plot3(iterations,0.1*oneD_arr,L2_C2_hp1,'--*m','LineWidth', 2);
h3 = plot3(iterations,0.05*oneD_arr,L2_C2_hp05,'--*b','LineWidth', 2);
h4 = plot3(iterations,0.025*oneD_arr,L2_C2_hp025,'--*c','LineWidth', 2);
h5 = plot3([40  40  40  40],h,1.2*h.^2,'-k','LineWidth', 2);
plot3([100 100 100 100],h,1.2*h.^2,'-k','LineWidth', 2);
plot3([160 160 160 160],h,1.2*h.^2,'-k','LineWidth', 2);
grid on;
set(gca, 'YScale', 'log', 'ZScale', 'log');
xlabel('decoupled iterations','Interpreter','latex'); ylabel('h','Interpreter','latex'); zlabel('$L_2$','Interpreter','latex');
legend([h1(1),h2(1),h3(1),h4(1),h5(1)], 'h = 0.2','h = 0.1','h = 0.05','h = 0.025','slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);
view(176,20.6);

% fig = figure(2);
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'iterations_error_tradeoff','-dpdf','-r0')

%% problem parameters
                   
% implicit hole cutting options -- apply to both grids
% ov_info = containers.Map({'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                           'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order', 'solve type', 'fringe update', 'exchange iterations'}, ...
%                          {2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "radial", ...
%                           2.5*max(h2), "rbf", "gaussian", 1.0, 2, "decoupled", "direct", 25});
% 
% % time step and linear solve parameters
% lin_sol_info = containers.Map({'type', 'Newton steps', 'Newton tolerance', 'GMRES tolerance', 'GMRES iterations', ...
%                                'decoupled loops', 'solution tolerance', 'dsol tolerance'}, ...
%                               { "iterative", 1, 1e-09, 1e-06, 100, 200, 1e-5, 1e-9});
                          
del    = 0.2;
h      = [del del/2 del/4 del/8];
L2_C2  = [5.152661e-02 9.918104e-03 2.935698e-03 7.070890e-04];
L2_C2_1em4_dir = [5.1542861272801152e-02,9.9165650735241716e-03,2.9457114181139889e-03,8.7448687582701817e-04];
L2_C2_1em5_dir = [5.1529189755299389e-02,9.9178520553963329e-03,2.9366559652752172e-03,7.2237728382130149e-04];

y_vals2= 0.1*h.^2;

figure(3)
clf
hold on
plot(h, L2_C2,'--*r','LineWidth', 2);
plot(h, L2_C2_1em5_dir,'--*b','LineWidth', 2);
plot(h,y_vals2,'-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('RBF$(C^2)$, coupled','RBF$(C^2)$, tol=1e-05','slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

%% problem parameters

% ov_info = containers.Map({'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                           'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order', 'solve type', 'fringe update', 'exchange iterations'}, ...
%                          {2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "radial", ...
%                           2.5*max(h2), "rbf", "gaussian", 1.0, 2, "decoupled", "iterative", 250});
% 
% % time step and linear solve parameters
% lin_sol_info = containers.Map({'type', 'Newton steps', 'Newton tolerance', 'GMRES iterations', 'GMRES tolerance', ...
%                                'decoupled loops', 'solution tolerance', 'dsol tolerance'}, ...
%                               { "iterative", 1, 1e-07, 1, 1e-15, 1, 1e-5, 1e-9});

L2_C2           = [5.1526623792468251e-02, 9.9181017536296504e-03, 2.9357283922972262e-03, 7.0705935537036422e-04];
L2_C2_hp2_itr   = [0.104936051851014,0.0675715470054935,0.0654152856793051,0.0561759221623225,0.0558805164393690,0.0530539526336566,0.0529654444372648,0.0520588243740360,0.0520134841895668,0.0517154012210150,0.0516928860922830,0.0515937676718178,0.0515836150529602,0.0515504618543996,0.0515461917440647,0.0515350602735441,0.0515333458089681,0.0515295974018441,0.0515289305504736,0.0515276652569849];
L2_C2_hp1_itr   = [0.0847062576328684,0.0522000988156351,0.0373599998265917,0.0265886051847968,0.0201512314212890,0.0151033668882459,0.0128417028468156,0.0107762874782609,0.0103445611699660,0.00964064793661243,0.00972445882710858,0.00950888737188921,0.00965009094364340,0.00959150687658908,0.00969735162452013,0.00968840615585061,0.00975682964442401,0.00976271174006895,0.00980524622378034,0.00981429947087820,0.00984060978239835,0.00984901976750969,0.00986541399275239,0.00987216178410443,0.00988249301996633,0.00988755049982216,0.00989413744247973,0.00989778482356542,0.00990202848192106,0.00990459581231119,0.00990735260204164,0.00990913024357180,0.00991093162200952,0.00991214821666230,0.00991332953270045,0.00991415512619401,0.00991493133850664,0.00991548811074876,0.00991599865954952,0.00991637241076241,0.00991670840850139,0.00991695843805324,0.00991717964260744,0.00991734647244314,0.00991749214580170,0.00991760324192319,0.00991769920081985,0.00991777307013393,0.00991783629768881,0.00991788535629044];
L2_C2_hp05_itr  = [0.0806980743896907,0.0556702066530226,0.0473999343777435,0.0380923018568279,0.0332446897608421,0.0275897016627544,0.0244051906242420,0.0207576684067476,0.0186060329587552,0.0161619082189809,0.0146676960540864,0.0129803175971140,0.0119146444375055,0.0107165913162262,0.00993615414674886,0.00906131525546557,0.00847473989215029,0.00781835369935319,0.00736691693807433,0.00686239828514595,0.00650812845193540,0.00611269921351799,0.00583066942714406,0.00551626481648436,0.00528961707228117,0.00503719972956602,0.00485404173677784,0.00465014340125502,0.00450169146181043,0.00433636751138179,0.00421586922900688,0.00408150767668472,0.00398362235974300,0.00387424537478186,0.00379467833640522,0.00370551726878207,0.00364078929020387,0.00356800891962811,0.00351529672468369,0.00345580197316666,0.00341281743143685,0.00336410818900355,0.00332900102673093,0.00328905783354958,0.00326033495995092,0.00322752674097726,0.00320398493110789,0.00317699389968951,0.00315766393644235,0.00313542450958341,0.00311952529625598,0.00310117462793430,0.00308807588089869,0.00307291401273819,0.00306210619746609,0.00304956409980800,0.00304063436452941,0.00303024838355030,0.00302286137773260,0.00301425283509244,0.00300813545567070,0.00300099436278937,0.00299592362961638,0.00298999567202691,0.00298578908143954,0.00298086519535719,0.00297737302980371,0.00297328103223469,0.00297038021224518,0.00296697805936112,0.00296456721248162,0.00296173755174500,0.00295973304363347,0.00295737879688399,0.00295571152537492,0.00295375229631458,0.00295236509267856,0.00295073423605387,0.00294957975140161,0.00294822197533629,0.00294726095624654,0.00294613035577535,0.00294533023347095,0.00294438867407941,0.00294372240718541,0.00294293819340570,0.00294238331655519,0.00294173009343410,0.00294126793354342,0.00294072377887809,0.00294033880841876,0.00293988548165368,0.00293956478417053,0.00293918710422359,0.00293891993204509,0.00293860526144099,0.00293838266927468,0.00293812048607459,0.00293793502718148,0.00293771656952346];
L2_C2_hp025_itr = [0.0772257299644825,0.0564430209736390,0.0490720449921352,0.0463622360372853,0.0448116933558921,0.0432942861872096,0.0419301737626522,0.0404024189615006,0.0389562063486881,0.0374143326505026,0.0359656558864569,0.0344832367773484,0.0331087396715319,0.0317321503372179,0.0304656354708379,0.0292081109932605,0.0280552175931101,0.0269134041593592,0.0258680700067359,0.0248325566359638,0.0238849155625577,0.0229448258831439,0.0220844278235785,0.0212292341665316,0.0204462728882321,0.0196664379758464,0.0189521525600506,0.0182392839284814,0.0175860192398029,0.0169328269288720,0.0163339590027772,0.0157341415419624,0.0151839570103559,0.0146320696608875,0.0141256377625481,0.0136169656526115,0.0131500188142324,0.0126804635585310,0.0122492922909285,0.0118152743801851,0.0114166375999324,0.0110150140921235,0.0106460639632598,0.0102740599091120,0.00993227952443502,0.00958743123192291,0.00927058326383555,0.00895069153014558,0.00865677605895344,0.00835986767094387,0.00808708758731941,0.00781138361481852,0.00755811621591630,0.00730200584952407,0.00706677946802871,0.00682879836332769,0.00661027306902290,0.00638908534407416,0.00618603693224535,0.00598042012317787,0.00579172675246859,0.00560055916418471,0.00542519013009896,0.00524744012828518,0.00508444636013111,0.00491916323999481,0.00476766919177258,0.00461397533129349,0.00447317200383617,0.00433025606767658,0.00419939495298504,0.00406650596003022,0.00394489374159008,0.00382133584654219,0.00370832972350597,0.00359345759282071,0.00348846112298966,0.00338167580922274,0.00328413518717295,0.00318488042470930,0.00309428112832639,0.00300203998849416,0.00291790374053975,0.00283219559605697,0.00275407759793392,0.00267445535614793,0.00260194175935682,0.00252798933131282,0.00246069491866967,0.00239202489709537,0.00232959095098394,0.00226584247511435,0.00220793481416158,0.00214877160319649,0.00209507877201313,0.00204018731209107,0.00199041891128338,0.00193950678331421,0.00189339192897833,0.00184618626661350,0.00180347217007016,0.00175971823859191,0.00172016889828808,0.00167962878633197,0.00164302378468387,0.00160547520153095,0.00157160860006759,0.00153684377180009,0.00150552309833063,0.00147334775648462,0.00144439307820333,0.00141462553472031,0.00138786861122416,0.00136033891354527,0.00133562242371830,0.00131017158383109,0.00128734842048780,0.00126382771167086,0.00124276033779898,0.00122103065276052,0.00120159051319256,0.00118152177730638,0.00116358875971838,0.00114505939315254,0.00112852133245584,0.00111141775518926,0.00109616997566036,0.00108038614969630,0.00106633103957882,0.00105176804309238,0.00103881465689401,0.00102538028556998,0.00101344396984340,0.00100105236125871,0.000990054400272530,0.000978625677813486,0.000968492956155935,0.000957952889610802,0.000948617569395605,0.000938897249992951,0.000930296460924692,0.000921331989171665,0.000913407530255337,0.000905139715447652,0.000897837767578095,0.000890211838289204,0.000883482687322082,0.000876448012527638,0.000870245782711059,0.000863755603464572,0.000858038001304769,0.000852049173056701,0.000846777241808656,0.000841249987564305,0.000836387872590294,0.000831285547140583,0.000826800272378966,0.000822089137828828,0.000817950393574457,0.000813599406345518,0.000809779348474576,0.000805759957886704,0.000802233018572601,0.000798518977021897,0.000795261686895913,0.000791828871553969,0.000788819693170056,0.000785645939037142,0.000782865111412455,0.000779930055472963,0.000777359449395795,0.000774644385547788,0.000772267369277310,0.000769755113644372,0.000767556428571640,0.000765231194754587,0.000763196840551995,0.000761044124338538,0.000759161253095223,0.000757167726120427,0.000755424544943879,0.000753577956778404,0.000751963638333457,0.000750252726470515,0.000748757326928737,0.000747171734144719,0.000745786118025440,0.000744316316597561,0.000743032087994697,0.000741669310274971,0.000740478749984625,0.000739214924848074,0.000738110932931788,0.000736938627640327,0.000735914670984739,0.000734827038031698,0.000733877102490459,0.000732867831017126,0.000731986377748150,0.000731049649149864,0.000730231574786957,0.000729362022153029,0.000728602622477696,0.000727795293534552,0.000727090230339214,0.000726340553592756,0.000725685824451559,0.000724989578246915,0.000724381488932613,0.000723734773174473,0.000723169912483550,0.000722569122780202,0.000722044339547502,0.000721486143562983,0.000720998525667066,0.000720479841484264,0.000720026696662864,0.000719544672977538,0.000719123511288270,0.000718675509271249,0.000718284027047858,0.000717867603725797,0.000717503668895058,0.000717116561915302,0.000716778200550692,0.000716418314206879,0.000716103698208288,0.000715769090615909,0.000715476526414000,0.000715165397734309,0.000714893315927983,0.000714603997548558,0.000714350943393253,0.000714081887981849,0.000713846512652612,0.000713596285009165,0.000713377337571476,0.000713144606208926,0.000712940926235078,0.000712724455627002,0.000712534966402012,0.000712333609840993,0.000712157312107992,0.000711970005122237,0.000711805971320697,0.000711631725602892,0.000711479094558454,0.000711316992364094,0.000711174964391499,0.000711024153367894];

figure(4)
clf
sgtitle('$L_2$ error evolution with intermediate fringe update, 1 decoupled loop','Interpreter','latex', 'FontSize', 18)

subplot(2,2,1)
hold on
plot(1:length(L2_C2_hp2_itr), L2_C2(1)*ones(length(L2_C2_hp2_itr),1),'--*k','LineWidth', 2);
plot(1:length(L2_C2_hp2_itr), L2_C2_hp2_itr,'--*b','LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('GMRES iterations','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
title('$h=0.2$','Interpreter','latex');
legend('RBF$(C^2)$, coupled','RBF$(C^2)$, decoupled','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

subplot(2,2,2)
hold on
plot(1:length(L2_C2_hp1_itr), L2_C2(2)*ones(length(L2_C2_hp1_itr),1),'--*k','LineWidth', 2);
plot(1:length(L2_C2_hp1_itr), L2_C2_hp1_itr,'--*b','LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('GMRES iterations','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
title('$h=0.1$','Interpreter','latex');
legend('RBF$(C^2)$, coupled','RBF$(C^2)$, decoupled','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

subplot(2,2,3)
hold on
plot(1:length(L2_C2_hp05_itr), L2_C2(3)*ones(length(L2_C2_hp05_itr),1),'--*k','LineWidth', 2);
plot(1:length(L2_C2_hp05_itr), L2_C2_hp05_itr,'--*b','LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('GMRES iterations','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
title('$h=0.05$','Interpreter','latex');
legend('RBF$(C^2)$, coupled','RBF$(C^2)$, decoupled','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

subplot(2,2,4)
hold on
plot(1:length(L2_C2_hp025_itr), L2_C2(4)*ones(length(L2_C2_hp025_itr),1),'--*k','LineWidth', 2);
plot(1:length(L2_C2_hp025_itr), L2_C2_hp025_itr,'--*b','LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('GMRES iterations','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
title('$h=0.025$','Interpreter','latex');
legend('RBF$(C^2)$, coupled','RBF$(C^2)$, decoupled','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);


del              = 0.2;
h                = [del del/2 del/4 del/8];
GMRES_itr_coup   = [4, 7, 12, 19];
GMRES_itr_decoup = [10, 20, 50, 150];

figure(5)
clf
hold on
plot(h,   GMRES_itr_coup,'--*k','LineWidth', 2);
plot(h, GMRES_itr_decoup,'--*r','LineWidth', 2);
set(gca, 'XScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('GMRES iterations','Interpreter','latex');
legend('RBF$(C^2)$, coupled','RBF$(C^2)$, decoupled with intermediate fringe update','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);