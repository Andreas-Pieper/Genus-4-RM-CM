AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");

SetDebugOnError(true);
SetVerbose("Reconstruction", true);
SetVerbose("Theta", true);

//intrinsic ModuliToBigPeriodMatrixNoam(H, points : fix_signs := false) -> AlgMatElt
function ModuliToBigPeriodMatrixNoam(H, points : fix_signs := false)
//{ Modified version of the ModiliToBigPeriodMatrix. }
    prec := Min([Precision(Parent(elt)) : elt in points]);
    CC := ComplexFieldExtra(prec);
    // should they be purely imaginary?
    //assert &and[Abs(Re(p)) lt CC`epscomp : p in points];
    OH := Integers(H);
	B := Basis(OH);
    g := Degree(H);
    betas := [[CC | Evaluate(B[i], pl : Precision := prec+10) : pl in InfinitePlaces(H)] : i in [1..g]];
    Pi1 := Transpose(Matrix(CC, betas));
    Pi2 := DiagonalMatrix(points)*(Transpose(Pi1)^-1);
    //bigres := HorizontalJoin(Pi2, Pi1);
    //return Pi1, Pi2;
    return HorizontalJoin(Pi2, Pi1);
end function;
//end intrinsic;

//intrinsic ModuliToSmallPeriodMatrixNoam(F, points) -> AlgMatElt
function ModuliToSmallPeriodMatrixNoam(H, points)
//{ Given a list of taus and a totally real field, compute a small period matrix. }
    return SmallPeriodMatrix(ModuliToBigPeriodMatrixNoam(H, points));
end function;
//end intrinsic;

label := "2.2.12.1-578.1-c";
prec := 1000;
R<x> := PolynomialRing(RationalsExtra(prec));
Hextra := NumberFieldExtra(x^4 - 10*x^3 + 20*x^2 + 25*x - 25 : prec := prec);
print Hextra;
CC<I> := ComplexField(prec);

taus := [0.6588698727219766677711228069136081096079333717347327274312523853468482695467751136192678203515927934887588356374203229526726476498695675357062584850387616280808551154745361843974236860545836301098997600508116739343759337123641220677864761869510726865344683850701152628381413642835743082488471708109603434761404125623340053239231030037107777367480313441153857927259007721959535113000090338899346543296535857334419159504146181993039558831395093266993639318684479549383426286341365982236879302671296211395701759702574343337461510704386699842079610112256387315484860786869504794652174354718884786005129320075981557437747805867730306736563399457735055659660281494546941571764490342960699657150044287662542588264763779531589308901156522020848803758431239707828249203345379092578920322331098393384432385261703982437552189955449712371378910578899835052814151091346955831674780905804678152759473089994992160849088448895059783341352642174201581757627741813757111103757124305691204192756674869939033125580734633p1000 + 0.2173290998890068233316657900695669469616656896065302505680681381549351701854413177347147114316240521979047246616632946266940770802393093552320446443774443293685146648066707485010096036982534934344115605638176612380338167263579423485173330180154124132319889033555820803298409309817560518474004160803916684957303595164471941712609085200605301645483996729791629723663740968916692908730159994007033559732307548527321841178042624564125540484721447527952833980338120934760577945470158870274636572362446223140420411612006449865031637018486020898942160638467476708208848204265074110025310273667197109613969345982910020620540652974572934901324039265317780979456969138687199126726953147300677971426034206801393742947745525075647646719931298826230340016066832918052691927633397158082433522388529126541092035518608196776987197805480392191637137935961046418129842344857565754829931494945053029793095004527364637366979395398630155803077475204230420812032614675925591749660804245195634523308843466606395657148975976p1000*I, 0.4843349247075070079718007752770024986978337609366033273986780343964778706654140336049042720239412338575178236712088603070342305637750122477122412741721149092803571774829831512879809326115770287443259266881317609675845037880896794159262314570790579803989213421363081063773728171144604220896667421525595220515678296011361666320719210298034845627125524848474421320529373558626779753648859330314870203834744027793415859857320444685271945152971160735365964188949459221546448563854116107833210616500774431675830484612443722384105512338551979422275358498732078185793047851127100216691039018155541569526244486781440213410897469404688965081766727301371862226591916758600104254026183050888952513939328073313871236628402087924257719249038791633656262580355596976736844080337323339125398385989479565283795007389979948251796489909443242771290265849310178202275942977673246354917015625714133706558725047918492897206054206609325568255427794283408918383264658521327776405117437795890621264077663197837113385763153125p1000 + 0.04418272878546617806886139816766197219177737520390100074163550533109755145723864246618749786169113305154933012425112527392980026268440782278919181247691549540461039756248089511940774411565484035766691023554327560328364137682394451316176648651662807405120419869207612575790777237024521326203458354634261293272280972565329946765155874312158246306650347641429003992943709594426408739777171126920015452944251921276285495139009576017072200877530107822920297293664840911149783956316074261294755945222679318143563507094205339161104318664433887608469820650896688406094815552880280906567652599000253712128683823815175120785646375099270990254040116655950957039513849363058560753828686674778431503027440523934069816844588369254320093932798873747966950640803702527324316538279081369854756507786796139610080054947890583182889544811245567900892533906305365691224474939057781058744402959956457914681184054065600429075332560127275926452252354584120444115849236184917146887844822176640731690109509592920728030648862580p1000*I, 1.210303960187158001111722187435203383733761892389472957300474578115821393507719276069405621215446956544696970493429821041152752057265335927998882008537025693058378384956439414562649140662424017500913695640630864006946422614046351382781148482392431197615003728352832490478794596904452927059982928854955831644608791102111678170410708031072400031646320157931354362471579266050746325123256208298584958892605662557388498844198544391255337658064918178481855427712389059139792218604459608176484644113596882503539960719038675174546622895836234893612803058954255818211504786175742077316755367454108374355411242283946134740747320741376328340643310233922609200229331556269476392288347972464237487833181841861474975772545063541225420769571872002371668388857754455751117236458591677738154494750530379359976041352279999616465760806773619444614528136255095144219777167838753905246014515560701742465726652070480958962871170955152812128614853324621882324490431193224632119114808585133445694644380020456400412560297251p1000 + 0.3532760908520349367570460132084648004543081374525384337114594067939397217568679038890510087633969151476146606056514097706049091634426849266500511555374529637992497465692910837017565708532663862394058309557440196781016794668787102245654595977494149510023781476294527421857038475963143965965979671167819608240965883429439961470760124276332062887848339939198749449016118483025900547743048447226014317192829942378985316912353636171835308943621992154974115572073905115395193328299819925747175835308872691758614930562795640795316437166794793967371372678647276952946145668709230439660900383569379570819572731351252215729644647856415964939225662940439737315980760535154165896391697937967923231120893663080630419288465820211936055503308897679351045832740212374547483681786617548765954581085728046321100211111257319379747397646445629542972918560390053193495546629850740609980636245560282661828365640314318752400153725085320642982527129918462010123758615308715628662082570245254570989633955773920670586011733859p1000*I, 2.965881021506261227940274479346171230140282615572557067605324471906024189053280981224220595900347337357684226266629617216870441087714816713800424339224825413466861224952036507618382851648317422071053289323315658423524590589100983795457122721593827688285686672971018999615242748663384942154562751663208686470453562829727272243761728847790576779341042416571505085133096255062189400485774207167744853195547470536972009142399762325529464242372460529082728635423737308279796009414972863108568308282176580834542548776451828566697035274884166156483834641967931412615029197333159238606303524925081356964703081403814697929859998739516509937283827459995918136128687550121708653538811939831011647350325527861857935896484604922437881077704591477787489872131921190811922015778724261390768096247476750526142062931481919669440962204083668108909491818495758828053033779786098478592125545390657110621608897007785588963276322579510337594559079603453183685562755156005067931398625929995602727336439634618360348469835595p1000 + 0.6400815074597503773238873287887534553728519132843168813236468861702194257060832022133779169225058743056420456291840141180074342139173464627926977949300837189309669479572383425165458118943374496556089517464811675673069246308200780341517897936559146243048589494520114690115510521809940760530926872340618940157151239899884233055351083711603153221989545625108744425818340295723662366141630235157147752727695348724499596778227416502342555162560589736133822538495545585415541278446326839296161878302593377928389578462923015992849415741993045001133930653875050825117574116379066087898115509977641143106479510363019678355638581409751394331120754343556185276237077219447193244276456339700412332343843513562703036958871672681928204098420021800184249714558191714163261219642858616537753231378974325681692039175059611405850818432922124663438355971981118101836475281128808627594498159476791993536083160290110062764162795053690308329721686978060225856022285654582986619978547767096803371739139843766943449415570650p1000*I];

/*
  Ps := [ KMatrixSpace(ComplexField(80), 4, 8) |
  [ ComplexField(80) | 0.31285376687641523795631632747074170837934306120126464364518153149292336883171743p80*I, 0.25582143333225903152284270699508996930413606860235410188347402777367428199918678p80*I, -0.41936261253235256914933245087119226303272406176223888978298137821917149023171620p80*I, 0.22600835109429192846919869977818902341127530770141073508273488726395314530199107p80*I, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000p80, 0.72240891759366497462773677570477043686059942525900608866746258412292234087084214p80, 0.10437492884377012642314994133913231914029024545324322653201396141766187805193970p80, -0.41836507468221022903462017416022321541360740677868274119865911463306073400298414p80, -0.19496714707366338279651620858417047043647141331986163997865760178918062179403593p80*I, -0.35100541561225384284814995698022029473730112597098433866880172213402658397205170p80*I, 0.95448605326563043337225240436796521988525267039494315797342246526374874253068779p80*I, -0.83398972968450686397326474079959059626241904931963907489878676349454713823576544p80*I, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000p80, 4.2775910824063350253722632242952295631394005747409939113325374158770776591291578p80, 3.6595570936564401771676763899295914454190913949352310491970887931718171963102554p80, 0.56426310843252568442085967106330886225267986736139415479231324651727934554627679p80, 0.032892419008140543216561863911136222425227351128895238171841847356979466184623166p80*I, 0.096413185993272504536527231487449147137642423823418402329358081670196523326902645p80*I, -0.15356837899827127039606005980188092768016771363715233366286722321411714699625118p80*I, 0.20678943092564577356043121941225465341947233476713527266317646267998413757508738p80*I, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000p80, 6.2868416845373235435123378637414850199242979335535547998802000271339184901993296p80, 7.9048756732872183917169246981071231376446071133593176620156486498391789530182320p80, 6.1672353680590310133431109375141257926552467694965444457428324564795584943447326p80, 0.045478521292260938164984297003559702670032434188030405915101509325693656414885816p80*I, -0.11593390920717150126647196141297519200358789677917531615634446897735606827455530p80*I, 0.13210945015896359018684724735235720174231100841997822834070311629113633294560098p80*I, -0.058523656949999602050107216809928812548103722690923023506406252118931023782461474p80*I, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000p80, -1.2868416845373235435123378637414850199242979335535547998802000271339184901993296p80, 0.33119230421257130469224897062415309779601124625220806225524859557134197261957282p80, 0.68686659819065353127064956558278856050568076992074414066351341163622289411197475p80 ],

  [ ComplexField(80) | 1.1544565818242058418774665467553729627784482101606543067468725868703358945276136p80*I, 0.94400249813451356802343937064306164375405293314001524741831708209017051361797249p80*I, -1.5474831357878901585475418180308065689020044775639740667229378252198926721766086p80*I, 0.83398972968450686397326474079959059626241904931963907489878676349454713823576544p80*I, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000p80, 0.72240891759366497462773677570477043686059942525900608866746258412292234087084214p80, 0.10437492884377012642314994133913231914029024545324322653201396141766187805193970p80, -0.41836507468221022903462017416022321541360740677868274119865911463306073400298414p80, -0.052835426935468593237734992353391713441046289195337825686885446206656381726947526p80*I, -0.095121261550429678631773904587272513472513868182108527565715933127190737441727611p80*I, 0.25866244075052321625826364846337480720110186134199331546522328357268794567425704p80*I, -0.22600835109429192846919869977818902341127530770141073508273488726395314530199108p80*I, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000p80, 4.2775910824063350253722632242952295631394005747409939113325374158770776591291578p80, 3.6595570936564401771676763899295914454190913949352310491970887931718171963102554p80, 0.56426310843252568442085967106330886225267986736139415479231324651727934554627679p80, 0.0093089121512221786666462930806121262674405470101989198240822416567329322051913114p80*I, 0.027285979739243127043099117138047226163277048011987929994156992751876473730737754p80*I, -0.043461520691035215963474403077429235902000159652790842178515640065656738401783444p80*I, 0.058523656949999602050107216809928812548103722690923023506406252118931023782461478p80*I, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000p80, 6.2868416845373235435123378637414850199242979335535547998802000271339184901993296p80, 7.9048756732872183917169246981071231376446071133593176620156486498391789530182320p80, 6.1672353680590310133431109375141257926552467694965444457428324564795584943447326p80, 0.16069531583443826103976016116348662655154557159449612683064906499734462298805521p80*I, -0.40964472077366728442570077012045125203572400085051955299315302328191464205966240p80*I, 0.46679991377866605028523359843488303257824929066425348432666216482583526572901094p80*I, -0.20678943092564577356043121941225465341947233476713527266317646267998413757508737p80*I, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000p80, -1.2868416845373235435123378637414850199242979335535547998802000271339184901993296p80, 0.33119230421257130469224897062415309779601124625220806225524859557134197261957282p80, 0.68686659819065353127064956558278856050568076992074414066351341163622289411197475p80 ]
  ];
Pi := Ps[1];
//Pi := Ps[2];
Pi1, Pi2 := SplitBigPeriodMatrix(Pi);
Pi_rev := HorizontalJoin(Pi2, Pi1);
Pi_sm := SmallPeriodMatrix(Pi);
print "Siegel reducing";
Pi_sm_red, Q := SiegelReduction(Pi_sm);
*/
/*
  Q0 := ChangeRing(Q,BaseRing(Pi));
  A := Submatrix(Q0, 1,1,4,4);
  B := Submatrix(Q0, 1,5,4,4);
  C := Submatrix(Q0, 5,1,4,4);
  D := Submatrix(Q0, 5,5,4,4);
  Q1 := BlockMatrix(2,2, [Transpose(D), Transpose(B), Transpose(C), Transpose(A)]);
  Pi := Pi*Q1;
  Pi1, Pi2 := SplitBigPeriodMatrix(Pi);
  tau := Pi1^-1*Pi2;
*/

Pi_sm := ModuliToSmallPeriodMatrixNoam(Hextra, taus);
Q_CC, C_CC := ReconstructCurveG4(Pi_sm);
invs_CC := InvariantsGenus4Curves(Q_CC, C_CC : normalize := false);