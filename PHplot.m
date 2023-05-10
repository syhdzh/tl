function PHplot(m,L,sitas,PT,HT)
% 压力和膜厚绘图
L=0.27;
load WYZ.mat
WY1=WYZ{10,1};WY2=WYZ{10,2};WY4=WYZ{10,4};WY3=WYZ{10,3};WY3=7*10^(-5)-WY3;
% WY3=WYZ{6,3};
% WY3=roundn(WYZ{6,3},-5);
% P1=PT(:,45);
% load 130um.mat
[n,m]=size(WY1);fai=0;
% P2=PT(:,45);
sitas=(fai+(0:2*pi/n:2*pi-2*pi/n))*180/pi/12;
Ls=(0:m-1)*(L/m); %轴向
[SITA,LL]=meshgrid(sitas,Ls);
% plot(sitas,P1,sitas,P2)
%
subplot(2,2,1)
mesh(SITA,LL,WY1')
set(gca,'xtick',0:5:30);xlim([0 30]);xlabel('周向角度(°)'); ylabel('轴向长度(m)'); zlabel('水膜变形量(m)');
subplot(2,2,2)
mesh(SITA,LL,WY2')
set(gca,'xtick',0:5:30);xlim([0 30]);xlabel('周向角度(°)'); ylabel('轴向长度(m)'); zlabel('水膜变形量(m)');
subplot(2,2,3)
mesh(SITA,LL,WY3')
kk=WY3(1,1);
set(gca,'xtick',0:5:30);xlim([0 30]);xlabel('周向角度(°)'); ylabel('轴向长度(m)'); zlabel('水膜变形量(m)');
subplot(2,2,4)
mesh(SITA,LL,WY4')
set(gca,'xtick',0:5:30);xlim([0 30]);xlabel('周向角度(°)'); ylabel('轴向长度(m)'); zlabel('水膜变形量(m)');
% % subplot(1,2,1);
PT=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,2.55057397107442e-05,2.85301411778415e-05,2.89297749351237e-05,2.89394104387504e-05,2.89337618077088e-05,2.89634159632105e-05,2.89609794131062e-05,2.89931270132112e-05,2.90120365411005e-05,2.90185598412760e-05,2.90155258405657e-05,2.89962140983012e-05,2.90179010094742e-05,2.90093327943874e-05,2.90376978956112e-05,2.90130078166355e-05,2.90356568133234e-05,2.90216656655365e-05,2.90435458626668e-05,2.90372339568613e-05,2.90407360260128e-05,2.90445313674395e-05,2.90456261822276e-05,2.90268891215713e-05,2.89966248090408e-05,2.90806612220332e-05,2.90508313288849e-05,2.90286727837608e-05,2.90429708852116e-05,2.90298784227164e-05,2.90484475682421e-05,2.90597314674148e-05,2.90532182856463e-05,2.90376360521981e-05,2.90141863212241e-05,2.90237883052319e-05,2.90055651263709e-05,2.90241695527549e-05,2.89900276430435e-05,2.90004462275744e-05,2.89767696199212e-05,2.89875719581329e-05,2.89693794786670e-05,2.89621205272140e-05,2.89442810849955e-05,2.89340241855818e-05,2.88744567828722e-05,2.85261619898850e-05,2.55538721169127e-05,0;0,4.89682007790197e-05,5.48626021028912e-05,5.56500354820560e-05,5.56665880307682e-05,5.56541420465971e-05,5.57148320882411e-05,5.57085336609688e-05,5.57746591469963e-05,5.58133122158194e-05,5.58259107310124e-05,5.58218414765079e-05,5.57804128967969e-05,5.58235900986955e-05,5.58084116852500e-05,5.58648228251884e-05,5.58148167677829e-05,5.58609543733432e-05,5.58315780788472e-05,5.58757750285967e-05,5.58639417291787e-05,5.58705446947402e-05,5.58797533818680e-05,5.58812222612574e-05,5.58431249217719e-05,5.57810451342926e-05,5.59524206463694e-05,5.58911445561861e-05,5.58457874563982e-05,5.58753344667974e-05,5.58477445778004e-05,5.58865136626817e-05,5.59098385787901e-05,5.58958938336201e-05,5.58665344926692e-05,5.58168026524726e-05,5.58354034151693e-05,5.58007451041588e-05,5.58373417534180e-05,5.57683751364175e-05,5.57896766347995e-05,5.57407034387145e-05,5.57626917761534e-05,5.57269406587073e-05,5.57118775954924e-05,5.56772980512150e-05,5.56564151607850e-05,5.55377251083281e-05,5.48546574747613e-05,4.90659999163623e-05,0;0,7.08356912209476e-05,7.93058312262254e-05,8.04592847869125e-05,8.04773738500738e-05,8.04552211141653e-05,8.05497567206124e-05,8.05375993294161e-05,8.06405143160011e-05,8.06998342410278e-05,8.07180037464295e-05,8.07146330357982e-05,8.06469293728189e-05,8.07137359050152e-05,8.06932400019360e-05,8.07775677835557e-05,8.07001308807495e-05,8.07726205127349e-05,8.07245747579860e-05,8.07930174402992e-05,8.07759826704446e-05,8.07853708652796e-05,8.08020871067037e-05,8.08028078253335e-05,8.07439800157141e-05,8.06459486811843e-05,8.09125895932761e-05,8.08172212328089e-05,8.07458753719772e-05,8.07929102593537e-05,8.07488883731284e-05,8.08103603150704e-05,8.08467435978079e-05,8.08241206062029e-05,8.07823263125392e-05,8.07024661696521e-05,8.07315430584916e-05,8.06815511517722e-05,8.07356801431089e-05,8.06297166363753e-05,8.06642332304805e-05,8.05864005967281e-05,8.06213974570335e-05,8.05682138632030e-05,8.05449706815197e-05,8.04946828666487e-05,8.04629045422560e-05,8.02859307609041e-05,7.92937308943711e-05,7.09869554508054e-05,0;0,9.15336013939832e-05,0.000102146832820662,0.000103636077067492,0.000103648626345851,0.000103611569212790,0.000103744282413352,0.000103724099139109,0.000103867365552325,0.000103948230382297,0.000103971427992339,0.000103970564485492,0.000103870756695986,0.000103965800201959,0.000103940923869957,0.000104052941816165,0.000103944272646465,0.000104048516443371,0.000103976129025250,0.000104072668058824,0.000104049889976991,0.000104062158957282,0.000104088964905845,0.000104087328737745,0.000104005683661416,0.000103864707844751,0.000104239388853948,0.000104106157691971,0.000104004047201920,0.000104072105202759,0.000104009594089157,0.000104096914218404,0.000104147495482489,0.000104114629254192,0.000104061609534968,0.000103946456321608,0.000103989685385031,0.000103925099980393,0.000103996140787212,0.000103849349481563,0.000103901797622008,0.000103789082920924,0.000103840785716486,0.000103769389248168,0.000103738090247945,0.000103673035106580,0.000103630111360168,0.000103396039184834,0.000102130232472433,9.17430879958826e-05,0;0,0.000111449207296487,0.000123634901253312,0.000125423747417936,0.000125420391051792,0.000125361062813604,0.000125537764023987,0.000125507075664932,0.000125694984784069,0.000125798048843241,0.000125825484314497,0.000125829216372648,0.000125689490303892,0.000125820075049150,0.000125791366348313,0.000125930519963591,0.000125785125383361,0.000125929162755677,0.000125824175013115,0.000125954584059868,0.000125924685303500,0.000125940286733995,0.000125980393805987,0.000125974248371869,0.000125866992800005,0.000125672953642834,0.000126173568989142,0.000125997458372661,0.000125857674589830,0.000125951743443942,0.000125868541828015,0.000125985609566690,0.000126051621188110,0.000126006484417357,0.000125943629870306,0.000125786516345942,0.000125850069738585,0.000125771352588902,0.000125858384045018,0.000125665388594808,0.000125743075156342,0.000125586899584278,0.000125661008466384,0.000125569857916819,0.000125531236943411,0.000125452502530850,0.000125398378874938,0.000125108618934214,0.000123613322181829,0.000111724543173289,0;0,0.000223993356255107,0.000241989383276270,0.000245550739882061,0.000245301559895052,0.000245014464832728,0.000245557690507034,0.000245414908624624,0.000245947436484742,0.000246206709712484,0.000246264706096782,0.000246320042150835,0.000245818683306266,0.000246253944391796,0.000246189306733313,0.000246546695070011,0.000246069208955039,0.000246580053012362,0.000246165834080891,0.000246600480619515,0.000246490228177145,0.000246537511024242,0.000246690878550652,0.000246642454167932,0.000246337018237723,0.000245665528511209,0.000247242076909114,0.000246688242504013,0.000246226344116681,0.000246557982811363,0.000246300521423739,0.000246655441690942,0.000246840372393097,0.000246704736599413,0.000246588277332506,0.000246070842717612,0.000246325001029388,0.000246141347164444,0.000246368297913051,0.000245781777554741,0.000246120948697251,0.000245581531321716,0.000245883120562497,0.000245624781521908,0.000245553804912288,0.000245383311598513,0.000245266409289554,0.000244620433656876,0.000241929276419065,0.000224808364085116,0;0,3.75363544963445,4.06557741710974,3.93003682409028,3.86052821689653,3.84901597316074,3.79541075635461,3.75974622876314,3.67939173705026,3.63743057001405,3.65114338178417,3.64758061666898,3.64781003043388,3.56883252241878,3.56664262788689,3.58553927195438,3.57759237120737,3.56327776762441,3.58412873379493,3.57382308537542,3.55256709043365,3.55101685363571,3.52775095531996,3.53889209091814,3.57072789331756,3.60561563162092,3.52877877926587,3.52031944864091,3.55644344021553,3.55313088829353,3.55823634466323,3.51566466984586,3.50252244733202,3.55346564185563,3.57986051607754,3.60110045290881,3.55436940722740,3.57901125185458,3.62496557757493,3.64518598045567,3.66218248583919,3.71669245558562,3.74046212143577,3.75295900634300,3.79501646537804,3.82630875024975,3.88184288274765,3.98025275331139,4.07304810943921,3.70748679580238,0;0,6.69086675696383,7.28001888922991,7.04692773552532,6.90723714215602,6.89428374629351,6.80242332414034,6.71267895201039,6.58774836360417,6.52078296304938,6.54108989850929,6.52538340036240,6.53041073168589,6.38376266670005,6.36551867883966,6.42070546162875,6.40687633840315,6.36532649390248,6.41570225835000,6.38564063383621,6.36037392503843,6.34557661821846,6.30414637570011,6.33051347881932,6.38281898536340,6.44786788387397,6.32388164143294,6.29988470064726,6.37265885049333,6.37153417883271,6.35824490341593,6.29797402165952,6.28266599635298,6.36899576530162,6.40637568638317,6.44884575467780,6.35869457001579,6.38716013356880,6.49062349895578,6.52587383588304,6.53798569026603,6.64952017603060,6.67887582921858,6.71479920960943,6.77725595669225,6.83324648018846,6.94254576563095,7.11783316005686,7.29216299294665,6.62668693841078,0;0,8.90706164436462,9.72734800450260,9.43660045392632,9.22791053892900,9.21701653426729,9.10019542739612,8.95205042847999,8.80768162686879,8.72741483444625,8.75043591069340,8.71794772107911,8.72820637811133,8.52827281633443,8.48691930943199,8.58712611327550,8.56486739560878,8.49420177366349,8.57440644832140,8.52361695724177,8.50232950913575,8.47104311835698,8.41545170577268,8.45798351385586,8.52001252504233,8.60653724179166,8.46430088824809,8.41877785420153,8.52327964727897,8.52827063299858,8.48619967593542,8.42466519662574,8.41371128397635,8.52400769271396,8.56144328320627,8.62196809392928,8.49588468970979,8.51533502207466,8.67975551023468,8.72144730026819,8.71965511872256,8.88243374252792,8.90935415163236,8.97032780167891,9.04179559001368,9.11539969256888,9.27322010088409,9.50380963429089,9.74190701214724,8.84780228217695,0;0,10.4997337440020,11.5015453855870,11.1905174737098,10.9172883593535,10.9070801328561,10.7766759076682,10.5756635264780,10.4277235061780,10.3413297464182,10.3658974775876,10.3159308711091,10.3285954371518,10.0903088031551,10.0254575260705,10.1715743153044,10.1358891646706,10.0421362764135,10.1463760072041,10.0795807068737,10.0639045816016,10.0184837705337,9.95238212931083,10.0091863988098,10.0719055139446,10.1693832195250,10.0337954928326,9.96269686645342,10.0904247343920,10.1042227676495,10.0325108402807,9.97899559716561,9.97513669435251,10.1017583839857,10.1330711204161,10.2060083277931,10.0532400799584,10.0586773041035,10.2803498015651,10.3187239363932,10.3036022636317,10.5061107360687,10.5298915296134,10.6115874785384,10.6878481143785,10.7728397622441,10.9708987290520,11.2367779189884,11.5167703770203,10.4612721227444,0;0,11.5537368546172,12.6899142497146,12.3911759640581,12.0606591416935,12.0472088520034,11.9130027742910,11.6705398553056,11.5278060117640,11.4390330561508,11.4659182534196,11.4009282217053,11.4113353523824,11.1487513774537,11.0652976166069,11.2519443424939,11.1974565896980,11.0912114266328,11.2100513820800,11.1349515872952,11.1229834802331,11.0686594102139,10.9957442053169,11.0628203593158,11.1193916426011,11.2172661453152,11.1073574893798,11.0092225551808,11.1498615096444,11.1742800967566,11.0779592509032,11.0363634949477,11.0394713445691,11.1776913928202,11.2006915004530,11.2790725323123,11.1090623541730,11.1017447640122,11.3713359832104,11.3973340502126,11.3753524668227,11.6028753520904,11.6270630280411,11.7221964769276,11.8030106285149,11.8945471520808,12.1222685398013,12.4057837225372,12.7046939652693,11.5457577741801,0;0,12.1466895975588,13.3723570374019,13.1102335433715,12.7341392001142,12.7122294585950,12.5827524186538,12.3130549730159,12.1791362373286,12.0892873347095,12.1205017688043,12.0450730593434,12.0481649450392,11.7735796827753,11.6800526368844,11.8971919009482,11.8198220752946,11.7134927058421,11.8358108516779,11.7610695339898,11.7494067397944,11.6922407330584,11.6163102595179,11.6883885231191,11.7343676227656,11.8235057894706,11.7513124581887,11.6277633112970,11.7706132397413,11.8068582625264,11.6939751532942,11.6643236159862,11.6721957850410,11.8193870918755,11.8348022133552,11.9115230423491,11.7328491532211,11.7184257024332,12.0224975637355,12.0292029487593,12.0096353743534,12.2462192424767,12.2762356004405,12.3766561467030,12.4633056654160,12.5579774001071,12.8033849875125,13.0898018628944,13.3858956064757,12.1716690852539,0;0,12.3408053939755,13.6174812116778,13.4081982132243,13.0006971027416,12.9657744511239,12.8482948839256,12.5658827585562,12.4411489249654,12.3500573833391,12.3880079996428,12.3076977052856,12.2993322997413,12.0230021417984,11.9294571560552,12.1642488524013,12.0625097732361,11.9677568384117,12.0826657679728,12.0161206445042,12.0015143455459,11.9468174890484,11.8716374246775,11.9431326671481,11.9763066740151,12.0500274698005,12.0204417415633,11.8761680504659,12.0114154356745,12.0601954093991,11.9394529116742,11.9194337399325,11.9286674840106,12.0833210150665,12.0935651142707,12.1625182045325,11.9824310005937,11.9686364338141,12.2914250568198,12.2751176588978,12.2671887749886,12.4975945809386,12.5386055187068,12.6369722957157,12.7303763815989,12.8259749681442,13.0768052808355,13.3543903021556,13.6295199979306,12.3966301503407,0;0,12.1947748924775,13.4878084651537,13.3375160781800,12.9166889853151,12.8652208736403,12.7657683030253,12.4839256628159,12.3668735252369,12.2736152293303,12.3204954589326,12.2408766325689,12.2187304372216,11.9487321371614,11.8650954336975,12.1034124502855,11.9791229089593,11.9052334723788,12.0032395490400,11.9509624159595,11.9311027248157,11.8826237976180,11.8116615360401,11.8774691690496,11.8976087672283,11.9523057883782,11.9634197850427,11.8059417014163,11.9254803728568,11.9870140187613,11.8663521496792,11.8523942064824,11.8590423695121,12.0200778865185,12.0282197754851,12.0850868720447,11.9092311470597,11.9040876763831,12.2289485902796,12.1896400917296,12.2006924048030,12.4115779645812,12.4673124359785,12.5578695907688,12.6576026867553,12.7527630202769,12.9974419788328,13.2574770006363,13.4981328944742,12.2740599818388,0;0,11.7542934623446,13.0364409240280,12.9440589067853,12.5283057859144,12.4592669317926,12.3822368130251,12.1119316921778,12.0005698913943,11.9042360760309,11.9612917442551,11.8871131598311,11.8512235006282,11.5937389776946,11.5279977096208,11.7561231986225,11.6146488822768,11.5671069661024,11.6413403933578,11.6067874881689,11.5809044739401,11.5403730850831,11.4765187547383,11.5327606824457,11.5411377484761,11.5764330363820,11.6207672658899,11.4595798242540,11.5577346267067,11.6315697585247,11.5171600815625,11.5055501871331,11.5059151769824,11.6718500495240,11.6806928004245,11.7234449762556,11.5560115156016,11.5658973017996,11.8769587948118,11.8184594899304,11.8523536058241,12.0335626821522,12.1052477690712,12.1844088637706,12.2881697984815,12.3819342441564,12.6105695861847,12.8469566858026,13.0452858761256,11.8481335384702,0;0,11.0662684897223,12.3135390899207,12.2700449654970,11.8792935129532,11.7936031160362,11.7416360100051,11.4912202719959,11.3838425775908,11.2839894648150,11.3511890855280,11.2858963703064,11.2386699680795,10.9981850087998,10.9554563258542,11.1612411336753,11.0110717409480,10.9912398870933,11.0377892410656,11.0218230634060,10.9905392413924,10.9578559065082,10.9031735816585,10.9476005595110,10.9465493002830,10.9650366596467,11.0310486810351,10.8766705276466,10.9503838854177,11.0352519993144,10.9312858590270,10.9187876624711,10.9099072582156,11.0784904741757,11.0900188923085,11.1189431094764,10.9627866234507,10.9914388250168,11.2747209004364,11.2041227877075,11.2608256661204,11.4057594322094,11.4920573270037,11.5582187695318,11.6620652193605,11.7535193575172,11.9585327759460,12.1675101515661,12.3209262994052,11.1635677611484,0;0,10.1680893544759,11.3626070954620,11.3559744664864,11.0068383954959,10.9083575711782,10.8818218835309,10.6569744321345,10.5529632057508,10.4501858468751,10.5258035442451,10.4712662867993,10.4171205129888,10.1969414254842,10.1785240254341,10.3527299622106,10.2045917541233,10.2095085625714,10.2279785197249,10.2288043175087,10.1940413559270,10.1676278949637,10.1231883206179,10.1554351185794,10.1477958893704,10.1543568977019,10.2286111771847,10.0910973086461,10.1401242704469,10.2336443547312,10.1423112405856,10.1268159953962,10.1071175512148,10.2748423315802,10.2898664000104,10.3072683793297,10.1643368973002,10.2118408224401,10.4565852183040,10.3832661091863,10.4585697018334,10.5648023959366,10.6616738931380,10.7150743518781,10.8138789124266,10.9018102895642,11.0784435127009,11.2580268936669,11.3689545703188,10.2596758067704,0;0,9.10244615558871,10.2273367380376,10.2422999858809,9.94941552449174,9.84352260224885,9.84058409046681,9.64474409665095,9.54472420167793,9.44078064916714,9.52154905980849,9.47815544437311,9.42270676184491,9.22521699656253,9.22880838129463,9.36569220148089,9.23107535125779,9.25443934499953,9.24749524361756,9.26143532244136,9.22574728249772,9.20332279251258,9.16917535187744,9.19056043559116,9.17926716732196,9.18007021067431,9.24931030521597,9.13701057083958,9.16350482037098,9.26215306704104,9.18415029728112,9.16481291846770,9.13430517711553,9.29655658598041,9.31473666242593,9.32424462420055,9.19580776369905,9.25880971838779,9.45805806647738,9.39219875686510,9.47869865069094,9.54753649098182,9.64909545992079,9.69109106351310,9.77951623279639,9.86229535308957,10.0086092507483,10.1583084872070,10.2326452723896,9.17993515476074,0;0,7.90635328448521,8.94799056034928,8.97172608507150,8.74268223690516,8.63677861454250,8.65314394240733,8.48792772229686,8.39392541965276,8.29198819649505,8.37318775993310,8.34010850333891,8.28918272736198,8.11630578282935,8.13622507348774,8.23426228601728,8.12369759118810,8.15678133533541,8.13005780562907,8.15202996120205,8.11823083954514,8.09748865988757,8.07273981021912,8.08591720868818,8.07357694441585,8.07483558577377,8.12826132970888,8.04634567591103,8.05458842585179,8.15359626873484,8.08856676162128,8.06594095763644,8.02654578021122,8.17764340009019,8.19767708134639,8.20338216600103,8.09045594761068,8.16237169289766,8.31370865396989,8.26453423963575,8.35254138773078,8.38895024004662,8.48801619464589,8.52065003830928,8.59403232697606,8.66956297055063,8.78625354184671,8.90664199695380,8.95263325976730,7.96565267380233,0;0,6.62511241435548,7.56814484394621,7.58933056151957,7.42694101145349,7.32836609826932,7.35792272351861,7.22367839914705,7.13864270915782,7.04293282165565,7.11914044586498,7.09481735975823,7.05330258363972,6.90655698289931,6.93519996583841,6.99680521020595,6.91792602119489,6.95154329062371,6.91246437976573,6.93715956120098,6.90760117245340,6.88699177006692,6.87006326105937,6.87826326773429,6.86693886507855,6.87362244796619,6.90450372677715,6.85424876133844,6.84974286361966,6.94351458580075,6.89069300527670,6.86638697426959,6.82162493864972,6.95559885658175,6.97566628509106,6.98114610022857,6.88458116991261,6.95710131357274,7.06240476299809,7.03633637835474,7.11591780804281,7.12730270790030,7.21681086948166,7.24204886494563,7.29742632477933,7.36340690501923,7.45311893444051,7.54575354292774,7.57197426786432,6.66434283147920,0;0,5.30198430607096,6.13157919169151,6.14496098509240,6.04353420975567,5.95965572773571,5.99487925556291,5.89097726760371,5.81832516854676,5.73383079082606,5.79955571929750,5.78218819516098,5.75308715173037,5.63372971413101,5.66283130633639,5.69413340108416,5.65001888700142,5.67604559683667,5.63313055371928,5.65564505322707,5.63199921532522,5.61111633398236,5.60005052981987,5.60635574548165,5.59740001423833,5.61209440909657,5.61897871168132,5.59725572908547,5.58610793996766,5.66864739938058,5.62718597462579,5.60352724200065,5.55827989516807,5.66949503958357,5.68771027437350,5.69522591985738,5.61587878633160,5.68026373365608,5.74552909624676,5.74459442970332,5.80716083966828,5.80262525648546,5.87644693862008,5.89593789569269,5.93260401560664,5.98682553315916,6.05341180476487,6.12064331058879,6.13482142348053,5.32357609535940,0;0,3.98975291428272,4.68817563878914,4.69185183572126,4.64116952843874,4.57706477517289,4.61046663627284,4.53543604038648,4.47812089587147,4.40959790386970,4.46047683266040,4.44846963708119,4.43224626833428,4.34093193695311,4.36369589346435,4.37319649416897,4.36121856256665,4.37467549826830,4.33582671322341,4.35273277590273,4.33562976455719,4.31538384715036,4.30826293749486,4.31475973640842,4.30882754392617,4.33101520443319,4.31772872002020,4.31769259691899,4.30543925925879,4.37143135310461,4.34064688746187,4.32002028020820,4.27954484314683,4.36396160170509,4.37884473203385,4.38883807939748,4.32733617087527,4.37664552107624,4.41069531031901,4.43158663310353,4.47210522838176,4.46063165154556,4.51494691154163,4.52968268481532,4.54947723209507,4.59038626541174,4.63793294663860,4.68298029563425,4.69058272653662,3.99745319707978,0;0,2.74159110021997,3.29114564939950,3.28847309204982,3.27267931240110,3.23105146891541,3.25635993605163,3.20768191827985,3.16745214008581,3.11851509994138,3.15219471911336,3.14437616179094,3.13890960242207,3.07538489098578,3.08787656462828,3.08518833489745,3.09691278768850,3.09699925551730,3.06841216654592,3.07807920060297,3.06729253181961,3.04953855758579,3.04481025482028,3.05213486534858,3.04920826045216,3.07520699228605,3.05002971927499,3.06225281369548,3.05312681736076,3.09897379424584,3.07817269507104,3.06249597325369,3.03151623826015,3.08760266039002,3.09834734956422,3.10946683545963,3.06599812819319,3.09655924498235,3.10982692265772,3.14398913457560,3.16211373174952,3.15136635113926,3.18516959492037,3.19571888810434,3.20261296255739,3.22976115437088,3.26191252496730,3.28872864947374,3.29286023973624,2.74113969695985,0;0,1.61876577007152,2.00077274734449,1.99594195447720,1.99889975640984,1.97823169340386,1.99227511572818,1.96603914075795,1.94260879837796,1.91427560059827,1.93146168408606,1.92705148826340,1.92819586322588,1.89062462067500,1.89317247793352,1.88692581821478,1.90927307019005,1.90014398613281,1.88455176735736,1.88763142739091,1.88213660485393,1.86905508999945,1.86582582697390,1.87298967535828,1.87245865497985,1.89619581688836,1.86981357021897,1.88449813142060,1.88038652599390,1.90565347155978,1.89383670968256,1.88401535585319,1.86524515471623,1.89511364777441,1.90166695477706,1.91135474168066,1.88513778580356,1.89807470413920,1.90052590602353,1.93561211088782,1.93657153569885,1.93100064781991,1.94677679978199,1.95342447173943,1.95275145883116,1.96729829232805,1.98679937319815,1.99997293399051,2.00166803531951,1.61498074808005,0;0,0.683111343832724,0.881106681151793,0.878451146605438,0.884919178005162,0.879657921266867,0.883549066424552,0.874107766776507,0.864842323890770,0.854407762519933,0.859239995033301,0.857547152258362,0.860309684840890,0.844765643899204,0.841977048419059,0.838112477650802,0.855976231442027,0.845959028874551,0.841684031615063,0.840897241883298,0.839178768087768,0.832412334075151,0.830526628007295,0.835367934211307,0.836027935241828,0.850639712406734,0.833608841127579,0.842834597376990,0.842824341682182,0.851225638827960,0.846646550982177,0.842400812108294,0.835380152948150,0.845203276199478,0.848083681026106,0.853667517475485,0.842538434509499,0.843877069559540,0.843266962378267,0.865973521969464,0.859820293740670,0.859603966904308,0.863106535685156,0.866221937663669,0.863515919099078,0.868415448822329,0.877262443514727,0.881598021661613,0.881430220517385,0.680116404821007,0;0,3.08786276821088e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,8.75044326012457e-07,2.96448294409214e-06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,2.23461231596872e-07,1.22368246136741e-06,3.62605672360111e-07,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,3.21335652410225e-07,2.48354754409062e-07,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
PT=PT*10^4.5;
mesh(SITA,LL,PT')
set(gca,'xtick',0:5:30);xlim([0 30]);xlabel('周向角度(°)'); ylabel('轴向长度(m)'); zlabel('水膜压力(Pa)');
% % % subplot(1,2,2);
% mesh(SITA,LL,HT')
% set(gca,'xtick',0:60:360);xlim([0 360]);zlim([-inf inf]); xlabel('周向角度(°)'); ylabel('轴向长度(m)'); zlabel('水膜厚度(m)');
% saveas(gcf, '膜厚与压力分布.jpg')
end

