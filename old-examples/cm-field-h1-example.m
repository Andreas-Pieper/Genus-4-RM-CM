AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
Attach("~/github/Genus-4-RM-CM/CM/findg4.m");

Attach("~/FlintWrapper.m");
Attach("~/github/Genus-4-RM-CM/flint-stuff/schottky-fast-theta.m");
load "~/github/Genus-4-RM-CM/CM/full_proc.m";

load "algebraize.m";
SetVerbose("CMExp",2);
SetVerbose("Reconstruction",true);
SetVerbose("Theta",true);
SetDebugOnError(true);

function OliveToIgusa(I)
    J1 := -15*I[1];
    J2 := 45/8*(3*I[1]^2 - 25/2*I[2]);
    J3 := 15/8*(9/2*I[1]^3 - 25*I[1]*I[2] - 375/4*I[3]);
    J4 := (J1*J3-J2^2)/4;
    J5 := 81/16*(-3*I[1]^5 + 7625/256*I[1]^3*I[2] + 13125/256*I[1]^2*I[3] - 9375/128*I[1]*I[2]^2 - 28125/128*I[2]*I[3] - 28125/256*I[4]);
    return [J1,J2,J3,J4,J5];
end function;



prec := 300;
//8.0.1265625.1
//coeffs := [1, -1, 0, 1, -1, 1, 0, -1, 1];
//8.0.13140625.1
coeffs := [1, 3, 5, 3, 4, -3, 5, -3, 1];
QQ := RationalsExtra(prec);
CC<I> := QQ`CC;
eps := CC`epscomp;
R<x> := PolynomialRing(QQ);
f := R!coeffs;
K<nu> := NumberFieldExtra(f);
_, K0incl, inv := IsCMField(K);
K0, incl := Explode(K0incl);
taus := FullEnumerationG4(f);
tau := taus[2];
//inv := FindCurve(tau);
// reducing precision
tau0 := tau;
tau := ChangeRing(tau, CC);
print "checking if Jacobian";
tau_red := SiegelReduction(tau);
sch := SchottkyModularForm(tau_red : prec := (prec div 3));
print sch;
assert Abs(sch) lt 10^(-prec/4);
eqns := ReconstructCurveG4(tau_red : flint := true);
Q,C := Explode(eqns);
invs, wts := InvariantsGenus4Curves(Q,C);
invs0 := invs;
//invs := WPSMultiply(wts, invs, invs[4]^(-1/wts[4]));
//invs := WPSMultiply(wts, invs, invs[9]^(-1/wts[9]));
//AlgRecG4(invs : normalized := true);

invs0 := ChangeUniverse(invs, ComplexFieldExtra(Floor(prec/3)));
n := Normalize(invs0, wts);
_, i0 := Min([n[i] ne 0 select wts[i] else 100 : i in [1..#n]]);
n := WPSMultiply(wts, invs0, (CC!invs0[i0])^(-1/wts[i0]));
n := [Abs(n[i]) lt 10^(-Precision(invs0[1])) select 0 else n[i] : i in [1..#n]];
igu := OliveToIgusa(n);
igu_alg := AlgebraizedInvariants(igu, RationalsExtra(prec));

f_rec := HyperellipticPolynomials(HyperellipticCurveFromIgusaInvariants(igu_alg));
AlgRecG4(n : normalized := true);


//invs := WPSMultiply(wts, invs, invs[4]^(-1/wts[4]));
invs0 := WPSMultiply(wts, invs0, invs[9]^(-1/wts[9]));

// computing class field and recognizing invariants
AttachSpec("~/github/hilbertmodularforms/ModFrmHilD/spec");
G, mp := TotallyPositiveUnits(K0);
U, phiU := UnitGroup(ZZK : GRH := true);
//for i in [1..NumberOfGenerators(U)] do print [ EmbedExtra(phiU(U.1) : iota := inf) : inf in InfinitePlacesExtra(K) ]; end for;
gensU := Generators(U);
gensV := [ (phiU(gen)*inv(phiU(gen))) @@ mp2 : gen in gensU ];
V := sub< G | gensV >;
Q, pQ := quo< G | V >;

// trying to recognize invariants
polys := [];
for i->el in invs do
  if el ne 0 then
    m := MinimalPolynomial(ComplexField(100)!el,4);
    Append(~polys, <i, m>);
  end if;
end for;

// rescale to try to remove sqrt(-2)
invs1 := invs;
polys1 := polys;
invs := WPSMultiply(wts, invs1, Sqrt(CC!-2)^(-1/5));

polys := [];
for i->el in invs do
  if el ne 0 then
    m := MinimalPolynomial(ComplexField(100)!el,4);
    Append(~polys, <i, m>);
  end if;
end for;

// recognize invariants
L<z> := NumberFieldExtra(x^2 - x - 1);
v := InfinitePlaces(L)[1];
invs_L := [];
for i := 1 to #invs do
  if i in [el[1] : el in polys] then
    for pair in polys do
      ind, pol := Explode(pair);
      if ind eq i then
        rs := [el[1] : el in Roots(pol,L)];
        for r in rs do
          if Abs(Evaluate(r,v) - invs[ind]) lt 10^-10 then
            Append(~invs_L, r);
          end if;
        end for;
      end if;
    end for;
  else
    Append(~invs_L, L!0);
  end if;
end for;

/* invs_L :=
 
[ ext<K0|Polynomial(K0, [-1, -1, 1])> where K0 is RationalsExtra(100) |
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | -19832929/1291467969, 271565/47832147 ],
[ RationalField() | -93241887983/92822968803906, 269895745/572981288913 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | -1/2, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | -1/4, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | -4352/35937, 35/1331 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | -2176/35937, 35/2662 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 19832929/2582935938, -271565/95664294 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ],
[ RationalField() | 0, 0 ]
];

S0<[t]> := PolynomialRing(Parent(invs_L[1]), 6);
S<X,Y> := PolynomialRing(S0, 2);

f := t[6]*(X^5*Y-Y^6); // up to a constant it is our reconstructed sextic
v := &+[t[i]*MonomialsOfDegree(S,4)[i] : i in [1..5]]; // we search for a binary quartic that might work

invs_rk3, wgt := InvariantsGenus4CurvesRank3(f,v);

J := Ideal([invs_L[i]-invs_rk3[i] : i in [1..60]]);
RadicalDecomposition(J);

/*
Groebner basis:
[
t[1],
t[2],
t[3],
t[4],
t[5] + 1/10469*(2835*$.1 + 10221)*t[6]^4,
t[6]^5 + 1/1331*(945*$.1 - 4352)
]
]
Thus we find v and f, but we can do better by rescaling both to get rid of the annoying t[6]
//*/
/*
f := t[6]^6/NormalForm(t[6]^5, J)*(X^5*Y-Y^6);
v := NormalForm(t[5], J)*Y^4;

v := S!Evaluate(v, [X/t[6],Y/t[6]]);
f := S!Evaluate(f, [X/t[6],Y/t[6]]);

L<nu> := Parent(invs_L[1]);
P3 := ProjectiveSpace(L, 3);
S<x,y,z,t> := CoordinateRing(P3);

Q := x*z-y^2;
E := 10469*t^3-(2835*nu+10221)*z^2*t+(945*nu+3407)*(x^2*y-z^3);

inv, wgt := InvariantsGenus4Curves(Q,E);
WPSEqual(wgt, inv, invs_L); // indeed isomorphic
*/


//we derive an equation of the genus 4 curve

L<nu> := NumberFieldExtra(Polynomial(RationalsExtra(), [1,-1,-1]));
P3 := ProjectiveSpace(L, 3);
S<x,y,z,t> := CoordinateRing(P3);

Q := x*z-y^2;
E := 10469*t^3-(2835*nu+10221)*z^2*t+(945*nu+3407)*(x^2*y-z^3);

C := Curve(P3, [Q,E]);
End := GeometricEndomorphismRepresentation(C);

for e in End do           
Polredabs(MinimalPolynomial(e[2]));
end for;

ends_ZZ := [ChangeRing(el[2], ZZ) : el in End];

K := ext<Rationals()| MinimalPolynomial(ends_ZZ[5])>;
AlgQ:= MatrixAlgebra<Rationals(), 8 | [ChangeRing(mat, Rationals()): mat in ends_ZZ]>;
_, P := Diagonalisation(AlgQ);
L:= BaseRing(Parent(P));
_, m := IsSubfield(K, L);
elts := [((P * ends_ZZ[i]*P^-1)[1,1]) @@ m : i in [1..8]];
O := sub<RingOfIntegers(K)|elts>;
O;




MinimalPolynomial(ends_ZZ[#ends_ZZ]);
m := $1;
Discriminant(m);
Factorization($1);
O := EquationOrder(m);
UnitGroup(O);
ClassGroup(O);
PicardGroup(O);
c := Conductor(O);
Norm(c);
IsPrincipal(c);


// compare
/*
T := RiemannSurface(R.1^5 + 1, 3 : Precision := 300);
tau_right := SmallPeriodMatrix(T);
*/

/*
csQ := [ ComplexField(300) | 0.354972361568494614621479312091980693126305414485215234259746037762036976241337687183141657823406056323329727451974242370008681803398538448830501501769323222009490937803825311348898194088622772968770936814482376820960678718338535806938912725088926779823869567539657688897471111522210175818285904987257p300 - 0.104039210568620570218646400541411085364409239642995195985722197493562695953893300825382534394442970479096277169468275807032587767972708734997672285218915117145879977788240479979726020481082288099940875487269113665277479150950773509151413342433037403940365626913577553687824424515492669510902848383708p300*I, 1.26683931362866398090429963710587103698865573577852554845127718113673304303577128817564861633150695572903178621277152090292056573508227009329636532591775714127648290767125624873224146042113976571758044295995581922789470187574884212023691954462848118437007097346449289634518657326851490075320270087118p300 - 0.258211143630106052283903112991932137418956741668986720944517945484761961241498342581491168234873406842437677947672273718969268972117451628518421495225092801161412029650306898657836049843091365686611676474593285766480510260814282791459222624684525939880568227922501604336461179463591126257784019845132p300*I, 0.600979930141308232006658095144931489925551336930252543530774515566368869015573623208742646971017570161007630500887627602491044298811236093735520738778571652922659261840587854133581266606518414171755440083561350919374335467112009205131508850295485183283198981293262002967608108931686008914481315942745p300 - 0.188873891247613698707916047613819013250972061269962781077668065885215513867926772748180332376507811657096355105799296148217345262397180581542807166798117210306473665008305947573066823436106110462971740401459036615322491245676762466343167961059219166008363128621542766608409022919532230446017353572833p300*I, 0.540275072116465883226745554384772489177621038988489829211718445490813334629645738700576475726978610489798063183560515063239033714146284699760289045834445499335977965673776430823229648504005457141856776642810041848137085629552657991361038940562362472976292159626434311169447951703430090854562894610852p300 - 0.218235498049141589484635080040025772851102700970945474686352170263023086083251379414381173140048183703468733222483860773493340231247229350720283536633257895657591843150292230938044828954528927961963746758997112501755652752631734521929161093881324107541802276274064400256182981570204102782437736654318p300*I, 0.646617252700309789059122378110416454573930984813030296122585938063881571946693313310942467729713725973548065079742816828555034199659325723963642451644302130732004108166859957423044942827496919427997163412155787092011890133547861768342024521475543537220613904058672908681259817185373007290665134793036p300 + 1.15061375490500730509905223179183085544354584451768275081035752494313226642906706974104058178138654566065747686672739152068252191922937909378714072198013042223477529132745927455660841545667507522215552364929305339797452211201123150315834580371673658543505909597290340625519549582242133445979033508425p300*I, 0.933452528988112657855600864428163273585605238451936963134203893506284472958004189891110398294553064452308159026850657129526966445091243661014108693015134783528061284903710950555488173436630093988365631841281559552412217305478486946462521190176196713237760954146475239954625806264699874646410402874732p300 + 0.414157668626383507055767645730144188045178414399175573959160658008513599768471534126198102992451833546106126045290289746127465548334697712860806520200107305392066687041985034356318958927283599579303225099542215332591151856257074125029526357028126867724714759650420937841920144928182039560899738021938p300*I, 0.813351086623696677355793186462564428235672245102598204128824161798328493124326478285392463914816250659781424239680408055211245485408814535439324983911697383399615549592560704445334968906429683557802719905659877089023346850343657991056719864906363710246475304946723440612846768910442740641281767587552p300 + 0.751788324902379227240584380230328495965043922280984173247693905537064411781301018534506575349028112200587331747527366202843151447536038136368189075541622990594721362366723666485918036863500845590544962968267157458005399315307756784147929320052359558504981211052475396336071685480147441176160291859225p300*I, 0.212948732997339043240573795220055639732752297173210478690683844531242512919479899311584018088797101149225171220281414239147868572833388821546142325147194000017259453378970892208730480234272839722081635388196117910315576462559967072750723271976536717953910068343042801742563775713628221796894755974721p300 - 0.0484311841435576807532250627070637344993802473851160796419904267619808345957182591108666692557778175125514723408663578245518053677286837132819380486648839179741405183598652463296923311297861382401401314815762204504720170888951447609261776400792246691912593609145356576478188744650749861994708799020454p300*I, 0.347423249118660581712073330432823925827503061306788042899782429736602974394296289068221126263926566690572777340523642947935113640049422231126120340125507204844377227733731397556545410406631452129365252110761075423093098334470678237004383259366091694975685994493822415303908982736342209861846253278599p300 - 0.00492949092675970624013456918181321312378676175838260300764808324948503432237467254581520899797190630273043434989962184868910297753513995484874328947684876511981468789919592263692470684467961696169724691691256065632638638196597390030760843619840929059357380602821460707370627719206307244853736782475205p300*I, 0.145218076911938447521921258450160351318744029537317990788501271485261870809351676593187315158196049721335832917156448653523498224292507317087350451713618181711306889574615454026163996935805191559010625176732988394390622377939370199109345851681031171447373890579621279791174932652195373490409780668519p300 + 0.0811167402308883626573095467905624553596089127096170892588686375188533188222117560879714318627941796877777424737354288361634146485904318827378289617114835444825197787445979012282509399269931142175902098747257066161375474913736915598720760035191655859034874304896264171156600724357297065482035771435163p300*I ];
monsQ := [ PolynomialRing(ComplexField(300), 4) | X^2, X*Y, X*Z, X*T, Y^2, Y*Z, Y*T, Z^2, Z*T, T^2 ];
csC := [ ComplexField(300) | 60.6305578210873480770462765220523688514034215603409438657743100154043278747302027078952032303444249969072772980423130926538329011011688911735418785381557190155999087402358420302107019853187844002063251992034296950582330881312889387471404584559325093448426189949303863858278884600134616856918905793721p300 + 25.8616522340654811328522839692882049789448089270477170797779530034957818626182614443370663159651503606215626874600852892303514464858713725653189992444827277806115577136378191371179857512244398323866390274931201512496116721304627355341949013052236301622067905104782966683312830276005075002371905250855p300*I, -56.5228479959497866777729648030274150499869093365505117162715475903178307609264246890740879381832336447217938504180277849443708228656705996374964526783377985849841731121034592778248571039134237786925119748813188942333210959309036254923975167595999797576143072654009997910574597119875922280555844254524p300 + 5.28280179256040194493398208835201647001525257501387077887737307612251154978728365276159305073165241400874742917495919636207232334647373456933871846006710714560114092997302199652387709992248981424554821750044934314815357305456134941010488258820135807295358153969870406927948933482791841021927594805993p300*I, 112.028046900486051986714396612934980745152066396295204300218203661771992211108411259957866730844136621074096874796543638427726159811915519603607326311623875563091317276824529513455287680364017286833930193397722106551142918821095402397273414774863965675360437217746127199988778483621584273217484990432p300 + 15.5342727801889200568708825316732768121679534911657485811250150424809101831867423601289509447091833607580926647379169987827899466377483592068702720542780172200873639285611638950073055750849614405240150868017652021962899896477927334265724926063674914134643504944015126251553427163662395949418374537273p300*I, -12.6031839782352425275438589561503565912228519635579915857096963296640379155321397806527060976256699219550725085410554804335199348960453486914558989603999739748226324299165960764762336867618217148703860152077016495580609672495855945332562464874835761099844013555345824112752153403555396728793580052658p300 + 1.15773030240231255049891298493191368576151391336479481945663643009549578915344887238541236308269737241553672074205972610677435072414100737283950416547964494525719223829785908601778064909633909826531771418495505620013225804356622031342561971646044188355407290442389006874565219112421701806257973927102p300*I, 64.2714985383985694044046694204234350141682668688869161133042771240737654600205999855958679356126725746009434515964515481078669981215149093480599623787854429733497536280616848918997417000061327280447861404282325823071362020076617282371092310386966014891644750484180078147662994631960673566961672248554p300 - 24.0339184772194958250150458294536878810758718710101873006128756605664713019102098256824771185840306342287482800866549649994405818022103322268355756956863844681213004786505773702732076481908350868895964047394521586956763055640395758743167325208900174394213012462320226765342715827427438922698818933934p300*I, -44.8175512549436773759886569236234292998483381363781826948706437925060599463075333545633104740420071208934716602530264925592966098415690347923409996272858436162048408501675262019704722061019872487847736715340841190346856672681316500253063075215682357925956330926700657837468854986805619390619298184353p300 + 19.5105287872159064677979791212804200471056275634364448397250723810421431426371360860772188115030389429184226669971146120952034650033862163809395932567104568594014757569983121968716732184029288252019731934767556961048399575855080202205143284789967381959286172515612965211666873063226249760296616707929p300*I, 47.8544971791657868119943832836692767630614731864667950221980036952138829078939508337949110935474996256457641142740385806199380928486304413658430488426177005736574099907726349597787425680286630587776072800231795728411670108806667334583824114754444123430014274000591805762090446444578162599935205769910p300 - 43.6725588750364842416947730945937381571583573532541958806574484797727400795024518718078528782971528714803048891222820898753492859023709340876537791602048302616714815218797631694790766901742767818312166421872038936749395617746571966567503958288714304415334790239824162621846369335870103616235686353467p300*I, 75.4714388341014535922494924239041822440162988410466039963002858350737972258470426095748978771812586653016336632475786213442382806932220358834577826037207759938225438109644849811690845105640969810847864630868109314521934509208568139135690750229443555322907749863842147945424045592165016457359543988971p300 - 1.47936505968412647490725850846633512229213409586465625213130682026447806176393738555233806175110396541222918053814307428978081804347999808874734701690315938710029494765708631222035092296445970652123409316713900777119431752923202971557556889018246429815032563822754895681790824036658539269057226121718p300*I, 5.83974748316233623633964345585308779264878573061642731541153448402006942247282474058630317388835415613963553905791460167168008646183535294192572989816418318148492472975489085089934491453049522797267458765898392034600871566471319323165740217497631799709143434196551922265530892575765566042676298454765p300 + 11.2300954918915892259918386190144196371066473475028199736178757825107493549452949309045100820818105712227780383125466860845436235602352033027087305049309948911277912857482399806532292691164815625546533817931255425429590093658175803011370327688093684318321493672433639474549034779733783110804324479505p300*I, 11.5422208986965600140217007605397012651101268210423679244873616657737923519170526152990117567451999318226461023515698963465897996595278269853871750010295767733929546354616083399524398348429218120305847188287107716256320508744042589704378356418483392014423207456874077857220008333624599277543023044103p300 - 21.2390519391383699373212475576926930180978274099371867347485420470463082370400675381269416207578279954037942257215457330467668338453250824705482871949327601375399978616046462789196434712158730818206113742682073242292140025031285006123614488987938115150630971435436989479313533806120829360420565310048p300*I, -26.8938728590091065080723430356391348788059110816820809171790443269397349710598694116211864058212302157080412645772057136238522194609335701114100256302615054546145460499853827427736631132824099989654543262251611238337089546186390695503766293038262134461868680574894083891456709130864026188250665535900p300 + 0.535656771755004528425502012764197003134210494698801676516384506243345809473153282319946655666759122380623623371859374354057167766608725447235285757141521081182550150566223733879026576197997010773361819935829920798839794236164414133625930253176944466057367484331866312671093995616999886303039446229088p300*I, 29.6886164306615726344470920501591421480948499380404975418863533222388793164150696805826076856020497049539922934445218250155565109551608326237681911428635847752964741665086168834132819860632489196995985659818059637282430069647133809869973782878260959029114835952397544027850869042277238087224538308116p300 - 46.5574429073011777624704270853702331748757487439738614046621003080584908010882704516495713135494116943972971411548689758050909819445346112834997129724686355213019507461895234195947908347772581541853467307417193431579866411988290379024546207747263794998535718182481219755176079776062129110133031820087p300*I, -31.5275457918828742326458971013854454111641869214156463175973032569997413062262209511082210076144194968642239870254412536895689084230668172450499358562446363409111047977523727274503313896523758223572279899747235933692151760805840589066912656184677208619242722073262859404393130979359978186615990385851p300 + 3.41633770368014142815327966561890557335669606252172060461115234804777030368027968043509239915564683652154299746109328594686524030853570658313797955441318424304289595305395155851052119869144410818448997068621990383153510109640921299875535624945973113598200417272000152738281781908194294817157127019708p300*I, 7.40523027346580772580387832342687861694011241031551916573261649256491266054865273234389910749490814946197843338380607601810737227116711616714145510481655836027365797527368863720981205892594155612700706937273054640235152562594939325506902028893077844621956711603746553730983303831691505359798920391945p300 + 4.93007365464170827034063820163018458295449092944083035057529544708512143964735338941756555653059300600872661647051377078515680220837565646866381192761345833076694959846679103184433680687962909897779107723396446265342497176144651141445264903460917393117067633683589294648613043791708858882933175221254p300*I, -9.17879813545111274794956948492113125723256306901099788723971044326804484053573579177586966456648506109276902460845377996586089334132879375439161876513158932297449656182500863715168139384727166907182936505869183466460669325165625517704064134027320358795764799612783582014987412890329164348571815378632p300 - 28.7360781066513703926302427740348242290772054897100265024473712430209219963755800143946584695192751822006065716482730897326174435367730977221852664306390270023366753389925707297155127832883765274545838680638703211955509418184904443550367842497172797608131294784958215820983927310759822401009616043517p300*I, -10.3844323650513208098367414845534467629232592279690712261629909253663458712952694317430039568314282620459990849298461076167418018449809019136244803993278312823512356483322483363517425511785154132894019468676677531988772911235232975380553330768997796742322794719171851098508239185926825297475477172439p300 + 2.19952422146000576740385830290051252892648725508270560979025398223928086885093412747586929571829550220420988607728127096127762652347687949799340806287012604622671364072425877299717668809405824134138157870207639460911057660616212033476116396672180955723180452648248226412123686229636555410624804575189p300*I, 19.6826844947130423171845768951098435330360007379114385662846950852809670076181633905761879590797676047643021400789917354093986650787092732610907973588109026148518696759101920691079999579670374548710767651360466796132882155069096868636270441813184827848043309596872715590037163033630050915808787206801p300 - 1.38622819650863341666836653644102995481486928666290236525797547004114857714733688911107939235961585743845715379095643088963807845187670562072368329404395173381636996053438477744812008018124743943748865672243437991445724975362425054056089625240179249391656998373860323718534618623057722730983099399116p300*I, 7.01336480542805664003441809247491058806011328521688789549055862189915035447258481117133967893059345710234574065408108623146440906805957300415454782582444288172745020891617981531829618628115371318033541967529765335800529398204543279473440409391601352876892370796082053589293488289565079741031982733827p300 + 10.7746612829277953583675016906567363693291126875179996358200397620356143661490981476734305186559590923060636365676728324335949557669719097524663338524664508202948487867055403220980788275215028209742688507331482306084974296661059222372332958247093857098183494800403113834366861152130865105080970441713p300*I, -4.89025692226439576264643660182182827194274735301878610518444659399814351034252362224934239508453584820443387999831072811451123147903062046309219399558616834245122285868783941306575045251954305143848124891329067928467473799795353444425467823273742295801195960187251741265306563869402627060559022285691p300 - 1.49864083201603481803020727167676250032824159365135115730615426621651879571150978949360567199584345725360981696194733865865952999755992790209724489872245494781347581769953707927573207709230812063086099419993080593996720874137154529405148822790579155603162455515181589147578341531617073230129672064774p300*I, -0.0188863027421601303142934265246277880613615345390004785538913656666686858053394396737719004448658249816764212961657822825874155137987399244217914116509956238259601674513090562016188593592401509369558163610609382729124281895149327004693555243564668461110056959058909651634360671730681621208855565988488p300 + 0.464670694777613222345912432065029679042387585999775646834229216896094990695418401383181991673881040198817257476077330963394192904133474956745986156462921169406595287808214100616890235039465982420219059820018116601545527420762113292375787639088730902584297688628666468602052298505590704482442505476941p300*I ];
monsC := [ PolynomialRing(ComplexField(300), 4) | X^3, X^2*Y, X^2*Z, X^2*T, X*Y^2, X*Y*Z, X*Y*T, X*Z^2, X*Z*T, X*T^2, Y^3, Y^2*Z, Y^2*T, Y*Z^2, Y*Z*T, Y*T^2, Z^3, Z^2*T, Z*T^2, T^3 ];
*/
