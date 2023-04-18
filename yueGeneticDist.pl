##Comparison.pl
my %hash;
my $r='';
open VCF,"veronica/Documents/School/McMaster/OtherSoftware/mt_all_noRoot.min4.fasta" or die "Can't open file '$file' for read. $!";
for $rec (<VCF>){
    chomp($rec);
    if($rec=~/>(.*)$/){
        $id=$1;
        print $id,"\n";
    }else{
        print $id,"\n";
        $hash{$id}=$rec;
        print $rec,"\n";
    }
}
close(VCF);
my @arr=();
@arr=("DRR015089",
"DRR015092",
"DRR015093",
"DRR015096",
"DRR015097",
"DRR015098",
"DRR015099",
"DRR015100",
"DRR015101",
"DRR015102",
"DRR015105",
"DRR015106",
"DRR015107",
"DRR015108",
"DRR015109",
"DRR015110",
"DRR015111",
"DRR015112",
"DRR017535",
"DRR017536",
"DRR017537",
"DRR017538",
"DRR017539",
"DRR017540",
"DRR017541",
"DRR017542",
"DRR017543",
"DRR017544",
"DRR017545",
"DRR017546",
"DRR017547",
"DRR017548",
"DRR017549",
"DRR017550",
"DRR017551",
"DRR017552",
"DRR017554",
"DRR017555",
"DRR017559",
"DRR017560",
"DRR017561",
"DRR022915",
"DRR022916",
"DRR022917",
"DRR022918",
"DRR022919",
"DRR022920",
"DRR022921",
"DRR022922",
"DRR022923",
"DRR022924",
"DRR022925",
"DRR022926",
"DRR022927",
"DRR032500",
"DRR032624",
"DRR146814",
"DRR237582",
"DRR237583",
"DRR289666",
"DRR289667",
"DRR289668",
"DRR289669",
"DRR289670",
"DRR290031",
"DRR290032",
"DRR290033",
"DRR290034",
"DRR290035",
"DRR290036",
"DRR290037",
"DRR290038",
"DRR290039",
"DRR290040",
"DRR290041",
"DRR290042",
"DRR290043",
"DRR290044",
"DRR290045",
"DRR290046",
"DRR290047",
"DRR290848",
"DRR290849",
"DRR290850",
"DRR290851",
"DRR290852",
"DRR290853",
"DRR290854",
"DRR290855",
"DRR295873",
"DRR295874",
"DRR295875",
"DRR377377",
"DRR377378",
"ERR232423",
"ERR232428",
"ERR2863830",
"ERR2864483",
"ERR2864484",
"ERR2864485",
"ERR3415715",
"ERR769498",
"ERR769499",
"ERR769500",
"ERR769501",
"ERR769502",
"ERR769503",
"ERR769504",
"ERR769505",
"ERR769506",
"ERR769507",
"ERR769508",
"ERR769509",
"ERR769510",
"ERR769511",
"ERR769512",
"ERR769513",
"ERR769514",
"ERR769515",
"ERR769516",
"ERR769517",
"ERR769518",
"ERR769519",
"ERR769520",
"ERR769521",
"ERR9791571",
"ERR9791572",
"ERR9791573",
"ERR9791574",
"ERR9791575",
"ERR9791576",
"ERR9791577",
"ERR9791578",
"ERR9791579",
"ERR9791580",
"ERR9791581",
"ERR9791582",
"ERR9791583",
"ERR9791584",
"ERR9791585",
"ERR9791586",
"ERR9791587",
"ERR9791588",
"ERR9791589",
"ERR9791590",
"ERR9791591",
"ERR9791592",
"ERR9791593",
"ERR9791594",
"ERR9791595",
"ERR9791596",
"ERR9791597",
"ERR9791598",
"ERR9791599",
"ERR9791600",
"ERR9791601",
"ERR9791602",
"ERR9791603",
"ERR9791604",
"ERR9791605",
"ERR9791606",
"ERR9791607",
"ERR9791608",
"ERR9791609",
"ERR9791610",
"ERR9791611",
"ERR9791612",
"ERR9791613",
"ERR9791614",
"ERR9791615",
"ERR9791616",
"ERR9791617",
"ERR9791618",
"ERR9791619",
"ERR9791620",
"ERR9791621",
"ERR9791622",
"ERR9791623",
"ERR9791624",
"ERR9791625",
"ERR9791626",
"ERR9791627",
"ERR9791628",
"ERR9791629",
"ERR9791630",
"ERR9791631",
"ERR9791632",
"ERR9791633",
"ERR9791634",
"ERR9791635",
"ERR9791636",
"ERR9791637",
"ERR9791638",
"ERR9791639",
"ERR9791640",
"ERR9791641",
"ERR9791642",
"ERR9791643",
"ERR9791644",
"ERR9791645",
"ERR9791646",
"ERR9791647",
"ERR9791648",
"ERR9791649",
"ERR9791650",
"ERR9791651",
"ERR9791652",
"ERR9791653",
"ERR9791654",
"ERR9791655",
"ERR9791656",
"ERR9791657",
"ERR9791658",
"ERR9791659",
"ERR9791660",
"ERR9791661",
"ERR9791662",
"ERR9791663",
"ERR9791664",
"ERR9791665",
"ERR9791666",
"ERR9791667",
"ERR9791668",
"ERR9791669",
"ERR9791670",
"ERR9791671",
"ERR9791672",
"ERR9791673",
"ERR9791674",
"ERR9791675",
"ERR9791676",
"ERR9791677",
"ERR9791678",
"ERR9791679",
"ERR9791680",
"ERR9791681",
"ERR9791682",
"ERR9791683",
"ERR9791684",
"ERR9791685",
"ERR9791686",
"ERR9791687",
"ERR9791688",
"ERR9791689",
"ERR9791690",
"ERR9791691",
"ERR9791692",
"ERR9791693",
"ERR9791694",
"ERR9791695",
"ERR9791696",
"ERR9791697",
"ERR9791698",
"ERR9791699",
"ERR9791700",
"ERR9791701",
"ERR9791702",
"ERR9791703",
"ERR9791704",
"ERR9791705",
"ERR9791706",
"ERR9791707",
"ERR9791708",
"ERR9791709",
"ERR9791710",
"ERR9791711",
"ERR9791712",
"ERR9791713",
"ERR9791714",
"ERR9791715",
"ERR9791716",
"ERR9791717",
"ERR9791718",
"ERR9791719",
"ERR9791720",
"ERR9791721",
"ERR9791722",
"ERR9791723",
"ERR9791724",
"ERR9791725",
"ERR9791726",
"ERR9791727",
"ERR9791728",
"ERR9791729",
"ERR9791730",
"ERR9791731",
"ERR9791732",
"ERR9791733",
"ERR9791734",
"ERR9791735",
"ERR9791736",
"ERR9791737",
"ERR9791738",
"ERR9791739",
"ERR9791740",
"ERR9791741",
"ERR9791743",
"ERR9791744",
"ERR9791745",
"ERR9791746",
"ERR9791747",
"ERR9791748",
"ERR9791749",
"ERR9791750",
"ERR9791751",
"ERR9791752",
"ERR9791753",
"ERR9791754",
"ERR9791755",
"ERR9791756",
"ERR9791757",
"ERR9791758",
"ERR9791759",
"ERR9791760",
"ERR9791761",
"ERR9791762",
"ERR9791763",
"ERR9791764",
"ERR9791765",
"ERR9791766",
"ERR9791767",
"ERR9791768",
"ERR9791769",
"ERR9791770",
"ERR9791771",
"ERR9791772",
"ERR9791773",
"ERR9791774",
"ERR9791775",
"ERR9791776",
"ERR9791777",
"ERR9791778",
"ERR9791779",
"ERR9791780",
"ERR9791781",
"ERR9791782",
"SRR068950",
"SRR068951",
"SRR068952",
"SRR10592629",
"SRR10592630",
"SRR10592631",
"SRR10592632",
"SRR10592633",
"SRR10714201",
"SRR10714202",
"SRR10714203",
"SRR10714204",
"SRR10714205",
"SRR10714206",
"SRR10714207",
"SRR10714208",
"SRR10714209",
"SRR10714210",
"SRR10714211",
"SRR10714212",
"SRR10714213",
"SRR10714214",
"SRR10714215",
"SRR10714216",
"SRR10714217",
"SRR10714218",
"SRR10714219",
"SRR10714220",
"SRR10714221",
"SRR10714222",
"SRR10714223",
"SRR10714224",
"SRR10714225",
"SRR10714226",
"SRR10714227",
"SRR10714228",
"SRR10714229",
"SRR10714230",
"SRR10714231",
"SRR10714232",
"SRR10714233",
"SRR10714234",
"SRR10714235",
"SRR10714236",
"SRR10714237",
"SRR10714238",
"SRR10714239",
"SRR10714240",
"SRR10714241",
"SRR10714242",
"SRR10714243",
"SRR10714244",
"SRR10714245",
"SRR10714246",
"SRR10714247",
"SRR10714248",
"SRR10714249",
"SRR10714250",
"SRR10714251",
"SRR10714252",
"SRR10714253",
"SRR10714254",
"SRR10714255",
"SRR10714256",
"SRR10714257",
"SRR10714258",
"SRR10714259",
"SRR10714260",
"SRR10714261",
"SRR10714262",
"SRR10714263",
"SRR10714264",
"SRR10815942",
"SRR10815943",
"SRR10815944",
"SRR10815945",
"SRR10815946",
"SRR10815947",
"SRR10815948",
"SRR10815949",
"SRR10815950",
"SRR10815951",
"SRR10815952",
"SRR10815953",
"SRR10815954",
"SRR10815955",
"SRR10815956",
"SRR10815957",
"SRR10815958",
"SRR10815959",
"SRR10815960",
"SRR10815961",
"SRR10815962",
"SRR10815963",
"SRR10815964",
"SRR10815965",
"SRR10815966",
"SRR10815967",
"SRR10815968",
"SRR10815969",
"SRR10815970",
"SRR10815971",
"SRR10815972",
"SRR10815973",
"SRR10815974",
"SRR10815975",
"SRR10997241",
"SRR10997243",
"SRR10997245",
"SRR10997250",
"SRR10997252",
"SRR10997258",
"SRR10997260",
"SRR10997263",
"SRR10997269",
"SRR10997272",
"SRR10997273",
"SRR10997274",
"SRR10997276",
"SRR10997277",
"SRR10997278",
"SRR10997281",
"SRR10997282",
"SRR10997283",
"SRR10997284",
"SRR10997285",
"SRR10997287",
"SRR10997288",
"SRR10997295",
"SRR10997297",
"SRR10997299",
"SRR10997300",
"SRR10997302",
"SRR10997305",
"SRR10997306",
"SRR10997309",
"SRR10997311",
"SRR10997313",
"SRR10997317",
"SRR10997322",
"SRR10997326",
"SRR10997327",
"SRR11425549",
"SRR11425550",
"SRR11785055",
"SRR11785056",
"SRR11785057",
"SRR11785058",
"SRR11785059",
"SRR11785060",
"SRR11785061",
"SRR11785062",
"SRR11785063",
"SRR11785064",
"SRR11785065",
"SRR11785066",
"SRR11785067",
"SRR11785068",
"SRR11785069",
"SRR11785070",
"SRR11785071",
"SRR11785072",
"SRR11785073",
"SRR11785074",
"SRR11785075",
"SRR11785076",
"SRR11785077",
"SRR11785078",
"SRR11785079",
"SRR11785080",
"SRR11785081",
"SRR11785082",
"SRR11785083",
"SRR11785084",
"SRR11785085",
"SRR11785086",
"SRR11785087",
"SRR11785088",
"SRR11785089",
"SRR11785090",
"SRR11785091",
"SRR11785092",
"SRR11785093",
"SRR11785094",
"SRR11785095",
"SRR11785096",
"SRR11785097",
"SRR11785098",
"SRR11785099",
"SRR11785100",
"SRR11785101",
"SRR11785102",
"SRR11785103",
"SRR11785104",
"SRR11785105",
"SRR11785106",
"SRR11785107",
"SRR11785108",
"SRR11785109",
"SRR11785110",
"SRR11785111",
"SRR11785112",
"SRR11785113",
"SRR11785114",
"SRR11785115",
"SRR11785116",
"SRR11785117",
"SRR11785118",
"SRR11785119",
"SRR11785120",
"SRR11785121",
"SRR11785122",
"SRR11785123",
"SRR11785124",
"SRR11785125",
"SRR11785126",
"SRR11785127",
"SRR11785128",
"SRR11785129",
"SRR11785130",
"SRR11785131",
"SRR11785132",
"SRR11785133",
"SRR11785134",
"SRR11785135",
"SRR11785136",
"SRR11785137",
"SRR11785138",
"SRR11785139",
"SRR11785140",
"SRR11785141",
"SRR11785142",
"SRR11785143",
"SRR11785144",
"SRR11785145",
"SRR11785146",
"SRR11785147",
"SRR11785148",
"SRR11785149",
"SRR11785150",
"SRR11785151",
"SRR11785152",
"SRR11785153",
"SRR11785154",
"SRR11785155",
"SRR11785156",
"SRR11785157",
"SRR11785158",
"SRR11785159",
"SRR11785160",
"SRR11785161",
"SRR11785162",
"SRR11785163",
"SRR11785164",
"SRR11785165",
"SRR11785166",
"SRR11785167",
"SRR11785168",
"SRR11785169",
"SRR11785170",
"SRR11785171",
"SRR11785172",
"SRR11785173",
"SRR11785174",
"SRR11785175",
"SRR11785176",
"SRR11785177",
"SRR11785178",
"SRR11785179",
"SRR11785180",
"SRR11785181",
"SRR11785182",
"SRR11785183",
"SRR11785184",
"SRR11785185",
"SRR11785186",
"SRR11785187",
"SRR11785188",
"SRR11785189",
"SRR11785190",
"SRR11785191",
"SRR11785192",
"SRR11785193",
"SRR11785194",
"SRR11785195",
"SRR11785196",
"SRR11785197",
"SRR11785198",
"SRR11785199",
"SRR11785200",
"SRR11785201",
"SRR11785202",
"SRR11785203",
"SRR11785204",
"SRR11785205",
"SRR11785206",
"SRR11785207",
"SRR11785208",
"SRR11785209",
"SRR11785210",
"SRR11785211",
"SRR11785212",
"SRR11785213",
"SRR11785214",
"SRR11785215",
"SRR11785216",
"SRR11785217",
"SRR11785218",
"SRR11785219",
"SRR11785220",
"SRR11785221",
"SRR11785222",
"SRR11785223",
"SRR11785224",
"SRR11785225",
"SRR11785226",
"SRR11785227",
"SRR11785228",
"SRR11785229",
"SRR11785230",
"SRR11785231",
"SRR11785232",
"SRR1187382",
"SRR11977809",
"SRR11977810",
"SRR11977811",
"SRR11977812",
"SRR11977813",
"SRR11977814",
"SRR11977815",
"SRR11977816",
"SRR11977817",
"SRR11977818",
"SRR11977819",
"SRR11977820",
"SRR11977821",
"SRR11977822",
"SRR11977823",
"SRR11977824",
"SRR11977825",
"SRR11977826",
"SRR11977827",
"SRR11977828",
"SRR11977829",
"SRR11977830",
"SRR11977831",
"SRR11977832",
"SRR11977833",
"SRR11977834",
"SRR11977835",
"SRR11977836",
"SRR11977837",
"SRR11977838",
"SRR11977839",
"SRR11977840",
"SRR11977841",
"SRR11977842",
"SRR11977843",
"SRR11977844",
"SRR11977845",
"SRR11977846",
"SRR11977847",
"SRR11977848",
"SRR11977849",
"SRR11977850",
"SRR11977851",
"SRR11977852",
"SRR11977853",
"SRR11977854",
"SRR11977855",
"SRR11977856",
"SRR11977857",
"SRR11977858",
"SRR11977859",
"SRR11977860",
"SRR11977861",
"SRR11977862",
"SRR11977863",
"SRR11977864",
"SRR11977865",
"SRR11977866",
"SRR12081224",
"SRR12190174",
"SRR12190175",
"SRR12190176",
"SRR12190177",
"SRR12190178",
"SRR12190179",
"SRR12190180",
"SRR12190181",
"SRR12190182",
"SRR12190183",
"SRR12763076",
"SRR12763077",
"SRR12763078",
"SRR12763079",
"SRR12763080",
"SRR12763081",
"SRR12763082",
"SRR12763083",
"SRR12763084",
"SRR12763085",
"SRR12763086",
"SRR12763087",
"SRR12763088",
"SRR12763089",
"SRR12763090",
"SRR12763091",
"SRR12763092",
"SRR12763093",
"SRR12763094",
"SRR12763095",
"SRR12763096",
"SRR12763097",
"SRR12763098",
"SRR12763099",
"SRR12763100",
"SRR12763101",
"SRR12763102",
"SRR12763103",
"SRR12763104",
"SRR12763105",
"SRR12763106",
"SRR12763107",
"SRR12763108",
"SRR12763109",
"SRR12763110",
"SRR12763111",
"SRR12763112",
"SRR12763113",
"SRR12763114",
"SRR12763115",
"SRR12763116",
"SRR12763117",
"SRR12763118",
"SRR12763119",
"SRR12763120",
"SRR12763121",
"SRR12763122",
"SRR12763123",
"SRR12763124",
"SRR12763125",
"SRR12763126",
"SRR12763127",
"SRR12763128",
"SRR12825367",
"SRR12849127",
"SRR12849128",
"SRR12849129",
"SRR12849130",
"SRR12849131",
"SRR12849132",
"SRR12849133",
"SRR12849134",
"SRR12894193",
"SRR12894194",
"SRR12894195",
"SRR12894196",
"SRR12894197",
"SRR12894198",
"SRR12894199",
"SRR12894200",
"SRR12894201",
"SRR12894202",
"SRR12894203",
"SRR12894204",
"SRR12949926",
"SRR12949927",
"SRR12949928",
"SRR12949929",
"SRR13127873",
"SRR13127874",
"SRR13127875",
"SRR13127876",
"SRR13127877",
"SRR13127878",
"SRR13127879",
"SRR13127880",
"SRR13127881",
"SRR13127882",
"SRR13127883",
"SRR13127884",
"SRR13579282",
"SRR13579283",
"SRR13579284",
"SRR13579285",
"SRR13579286",
"SRR13579287",
"SRR13579288",
"SRR13579289",
"SRR13579290",
"SRR13579291",
"SRR13579292",
"SRR13579293",
"SRR13579294",
"SRR13579295",
"SRR13579296",
"SRR13579297",
"SRR13579298",
"SRR13579299",
"SRR13579300",
"SRR13579301",
"SRR13579302",
"SRR13579303",
"SRR13579304",
"SRR13579305",
"SRR13579306",
"SRR13579307",
"SRR13579308",
"SRR13579309",
"SRR13579310",
"SRR13579311",
"SRR13579312",
"SRR13579313",
"SRR13579314",
"SRR13579315",
"SRR13579316",
"SRR13579317",
"SRR13579318",
"SRR13579319",
"SRR13579320",
"SRR13579321",
"SRR13579322",
"SRR13579323",
"SRR13579324",
"SRR13579325",
"SRR13579326",
"SRR13579327",
"SRR13579328",
"SRR13579329",
"SRR13579330",
"SRR13579331",
"SRR13579332",
"SRR13579333",
"SRR13579334",
"SRR13579335",
"SRR13579336",
"SRR13579337",
"SRR13579338",
"SRR13579339",
"SRR13579340",
"SRR13579341",
"SRR13579342",
"SRR13579343",
"SRR13579344",
"SRR13579345",
"SRR13579346",
"SRR13579347",
"SRR13579348",
"SRR13579349",
"SRR13579350",
"SRR13579351",
"SRR13579352",
"SRR13579353",
"SRR13579354",
"SRR13579355",
"SRR13579356",
"SRR13579357",
"SRR13579358",
"SRR13579359",
"SRR13579360",
"SRR13579361",
"SRR13579362",
"SRR13579363",
"SRR13579364",
"SRR13579365",
"SRR13579366",
"SRR13579367",
"SRR13579368",
"SRR13579369",
"SRR13579370",
"SRR13579371",
"SRR13579372",
"SRR13579373",
"SRR13579374",
"SRR13579375",
"SRR13579376",
"SRR13579377",
"SRR13579378",
"SRR13579379",
"SRR13579380",
"SRR13579381",
"SRR13579382",
"SRR13579383",
"SRR13579384",
"SRR13579385",
"SRR13579386",
"SRR13579387",
"SRR13579388",
"SRR13579389",
"SRR13579390",
"SRR13579391",
"SRR13579392",
"SRR13579393",
"SRR13579394",
"SRR13579395",
"SRR13579396",
"SRR13579397",
"SRR13579398",
"SRR13579399",
"SRR13579400",
"SRR13579401",
"SRR13579402",
"SRR13579403",
"SRR13579404",
"SRR13579405",
"SRR13579406",
"SRR13579407",
"SRR13579408",
"SRR13579409",
"SRR13579410",
"SRR13579411",
"SRR13579412",
"SRR13579413",
"SRR13579414",
"SRR13579415",
"SRR13579416",
"SRR13579417",
"SRR13579418",
"SRR13579419",
"SRR13579420",
"SRR13579421",
"SRR13579422",
"SRR13579423",
"SRR13579424",
"SRR13579425",
"SRR13579426",
"SRR13579427",
"SRR13579428",
"SRR13579429",
"SRR13579430",
"SRR13579431",
"SRR13579432",
"SRR13579433",
"SRR13579434",
"SRR13579435",
"SRR13579436",
"SRR13579437",
"SRR13579438",
"SRR13579439",
"SRR13579440",
"SRR13579441",
"SRR13579442",
"SRR13579443",
"SRR13579444",
"SRR13579445",
"SRR13579446",
"SRR13579447",
"SRR13579448",
"SRR13579449",
"SRR13579450",
"SRR13579451",
"SRR13579452",
"SRR13579453",
"SRR13579454",
"SRR13579455",
"SRR13579456",
"SRR13579457",
"SRR13579458",
"SRR13579459",
"SRR13579460",
"SRR13579461",
"SRR13579462",
"SRR13579463",
"SRR13579464",
"SRR13579465",
"SRR13579466",
"SRR13579467",
"SRR13579468",
"SRR13579469",
"SRR14364143",
"SRR14364144",
"SRR14364145",
"SRR14364146",
"SRR14364147",
"SRR14364148",
"SRR14364149",
"SRR14364150",
"SRR14364151",
"SRR14364152",
"SRR14364153",
"SRR14364154",
"SRR14364155",
"SRR14364156",
"SRR14364157",
"SRR14364158",
"SRR14364159",
"SRR14364160",
"SRR14364161",
"SRR14364162",
"SRR14364163",
"SRR14364164",
"SRR14364165",
"SRR14364166",
"SRR14584215",
"SRR14584216",
"SRR14584217",
"SRR14584218",
"SRR14584219",
"SRR14584220",
"SRR14584221",
"SRR14584222",
"SRR14584223",
"SRR14584224",
"SRR14584225",
"SRR14584226",
"SRR14584227",
"SRR14584228",
"SRR14584229",
"SRR14584230",
"SRR14584231",
"SRR14584232",
"SRR14584233",
"SRR14584234",
"SRR14584235",
"SRR14584236",
"SRR14584237",
"SRR14584238",
"SRR14584239",
"SRR14584240",
"SRR14584241",
"SRR14584242",
"SRR14584243",
"SRR14584244",
"SRR14584245",
"SRR14584246",
"SRR14584247",
"SRR14584248",
"SRR14584249",
"SRR14584250",
"SRR14584251",
"SRR14584252",
"SRR14584253",
"SRR14584254",
"SRR14584255",
"SRR14584256",
"SRR14584257",
"SRR14584258",
"SRR14584259",
"SRR14584261",
"SRR14584262",
"SRR14584263",
"SRR14584264",
"SRR14584265",
"SRR14584266",
"SRR14584268",
"SRR14584269",
"SRR14584270",
"SRR14584271",
"SRR14584272",
"SRR14584273",
"SRR14584274",
"SRR14584275",
"SRR14584276",
"SRR15010315",
"SRR15010316",
"SRR15010317",
"SRR15010318",
"SRR15010319",
"SRR15010320",
"SRR15010321",
"SRR15010322",
"SRR15010324",
"SRR15010326",
"SRR15010327",
"SRR15010328",
"SRR15010329",
"SRR15010330",
"SRR15010331",
"SRR15010332",
"SRR15010333",
"SRR15010334",
"SRR15010335",
"SRR15010336",
"SRR15010337",
"SRR15010338",
"SRR15010339",
"SRR15010340",
"SRR15010341",
"SRR15010342",
"SRR15010343",
"SRR15010344",
"SRR15010345",
"SRR15010346",
"SRR15010347",
"SRR15010348",
"SRR15010349",
"SRR15010350",
"SRR15010351",
"SRR15010352",
"SRR15010353",
"SRR15010354",
"SRR15010355",
"SRR15010356",
"SRR15010357",
"SRR15010358",
"SRR15010359",
"SRR15010360",
"SRR15010361",
"SRR15010362",
"SRR15010363",
"SRR15010364",
"SRR15010365",
"SRR15010366",
"SRR15010367",
"SRR15010368",
"SRR15010369",
"SRR15010370",
"SRR15010371",
"SRR15010372",
"SRR15010373",
"SRR15010374",
"SRR15010375",
"SRR15010376",
"SRR15010377",
"SRR15010378",
"SRR15010379",
"SRR15010380",
"SRR15010381",
"SRR15010382",
"SRR15010383",
"SRR15010384",
"SRR15010385",
"SRR15010386",
"SRR15010387",
"SRR15010389",
"SRR15010390",
"SRR15010391",
"SRR15010392",
"SRR15010393",
"SRR15010394",
"SRR15010395",
"SRR15010396",
"SRR15010397",
"SRR15010398",
"SRR15010399",
"SRR15010400",
"SRR15010401",
"SRR15010403",
"SRR15010404",
"SRR15010405",
"SRR15010406",
"SRR15010407",
"SRR15010408",
"SRR15010409",
"SRR15010410",
"SRR1523528",
"SRR15597331",
"SRR15597332",
"SRR15597333",
"SRR159251",
"SRR159252",
"SRR159253",
"SRR159254",
"SRR16079626",
"SRR16079627",
"SRR16079628",
"SRR1614167",
"SRR1614168",
"SRR1614169",
"SRR1614170",
"SRR1614172",
"SRR16287627",
"SRR16287628",
"SRR16539748",
"SRR16539749",
"SRR16539750",
"SRR16539751",
"SRR16539752",
"SRR16539753",
"SRR16539754",
"SRR16539755",
"SRR16539756",
"SRR16539757",
"SRR16539758",
"SRR16539759",
"SRR16539760",
"SRR16539761",
"SRR16539762",
"SRR16539763",
"SRR16539764",
"SRR16539765",
"SRR16539766",
"SRR16539767",
"SRR16539768",
"SRR16539769",
"SRR16539770",
"SRR16539771",
"SRR16539772",
"SRR16539773",
"SRR16539774",
"SRR16539775",
"SRR16539776",
"SRR16539777",
"SRR16539778",
"SRR16539779",
"SRR16539780",
"SRR16539781",
"SRR16539782",
"SRR16539783",
"SRR16539784",
"SRR16539785",
"SRR16539786",
"SRR16539787",
"SRR16539788",
"SRR16539789",
"SRR16539790",
"SRR16539791",
"SRR16539792",
"SRR16539793",
"SRR16539794",
"SRR16539795",
"SRR16539796",
"SRR16539797",
"SRR16944265",
"SRR16944266",
"SRR16944267",
"SRR16944268",
"SRR16944269",
"SRR16944270",
"SRR16944271",
"SRR16944272",
"SRR16944273",
"SRR16944274",
"SRR16944275",
"SRR16944276",
"SRR16944277",
"SRR16944278",
"SRR16944279",
"SRR16944280",
"SRR16944281",
"SRR16944282",
"SRR16944283",
"SRR16944284",
"SRR16944285",
"SRR16944286",
"SRR16944287",
"SRR16944288",
"SRR16944289",
"SRR16944290",
"SRR16944291",
"SRR16944292",
"SRR16944293",
"SRR16944294",
"SRR16944295",
"SRR16944296",
"SRR16944297",
"SRR16944298",
"SRR16944299",
"SRR16944300",
"SRR16944301",
"SRR16944302",
"SRR16944303",
"SRR16944304",
"SRR16944305",
"SRR16944306",
"SRR16944307",
"SRR16944308",
"SRR16944309",
"SRR16944310",
"SRR16944311",
"SRR16944312",
"SRR16944313",
"SRR16944314",
"SRR16944315",
"SRR16944316",
"SRR16944317",
"SRR16944318",
"SRR16944319",
"SRR16944320",
"SRR16944321",
"SRR16944322",
"SRR16944323",
"SRR16944324",
"SRR16944325",
"SRR16944326",
"SRR16944327",
"SRR16944328",
"SRR16944329",
"SRR16944330",
"SRR16944331",
"SRR16944332",
"SRR16944333",
"SRR16944334",
"SRR16944335",
"SRR16944336",
"SRR16944337",
"SRR16944338",
"SRR16944339",
"SRR16944340",
"SRR16944341",
"SRR16944342",
"SRR16944343",
"SRR16944344",
"SRR16944345",
"SRR16944346",
"SRR16944347",
"SRR16944348",
"SRR16944349",
"SRR16944350",
"SRR16944351",
"SRR16944352",
"SRR16944353",
"SRR16944354",
"SRR16944355",
"SRR16944356",
"SRR16944357",
"SRR16944358",
"SRR16944359",
"SRR16944360",
"SRR16944361",
"SRR16944362",
"SRR16944363",
"SRR16944364",
"SRR16944365",
"SRR16944366",
"SRR16944367",
"SRR16944368",
"SRR16944369",
"SRR16944370",
"SRR16944371",
"SRR16944372",
"SRR16944373",
"SRR16944374",
"SRR16944375",
"SRR16944376",
"SRR16944377",
"SRR16944378",
"SRR16944379",
"SRR16944380",
"SRR16944381",
"SRR16944382",
"SRR16944383",
"SRR16944384",
"SRR16944385",
"SRR16944386",
"SRR16944387",
"SRR16944388",
"SRR16944389",
"SRR16944390",
"SRR16944391",
"SRR16944392",
"SRR16944393",
"SRR16944394",
"SRR16944395",
"SRR16944396",
"SRR16944397",
"SRR16944398",
"SRR16944399",
"SRR16944400",
"SRR16944401",
"SRR16944402",
"SRR16944403",
"SRR16944404",
"SRR16944405",
"SRR16944406",
"SRR16944407",
"SRR16944408",
"SRR16944409",
"SRR16944410",
"SRR16944411",
"SRR16944412",
"SRR16944413",
"SRR16944414",
"SRR16944415",
"SRR16944416",
"SRR16944417",
"SRR16944418",
"SRR16944419",
"SRR16944420",
"SRR16944421",
"SRR16944422",
"SRR16944423",
"SRR16944424",
"SRR16944425",
"SRR16944426",
"SRR16944427",
"SRR16944428",
"SRR16944429",
"SRR16944430",
"SRR16944431",
"SRR16944432",
"SRR16944433",
"SRR16944434",
"SRR16944435",
"SRR16944436",
"SRR16944437",
"SRR16944438",
"SRR16944439",
"SRR16944440",
"SRR16944441",
"SRR16944442",
"SRR16944443",
"SRR16944444",
"SRR16944445",
"SRR16944446",
"SRR16944447",
"SRR16944448",
"SRR16944449",
"SRR16944450",
"SRR16944451",
"SRR16944452",
"SRR16944453",
"SRR16944454",
"SRR16944455",
"SRR16944456",
"SRR16944457",
"SRR16944458",
"SRR16944459",
"SRR16944460",
"SRR16944461",
"SRR16944462",
"SRR16944463",
"SRR16944464",
"SRR16944465",
"SRR16944466",
"SRR16944467",
"SRR16944468",
"SRR16944469",
"SRR16944470",
"SRR16944471",
"SRR16944472",
"SRR16944473",
"SRR16944474",
"SRR16944475",
"SRR16944476",
"SRR16944477",
"SRR16944478",
"SRR16944479",
"SRR16944480",
"SRR16944481",
"SRR16944482",
"SRR16944483",
"SRR16944484",
"SRR16944485",
"SRR16944486",
"SRR16944487",
"SRR16944488",
"SRR16944489",
"SRR16944490",
"SRR16944491",
"SRR16944492",
"SRR16944493",
"SRR16944494",
"SRR17181022",
"SRR17181023",
"SRR17181024",
"SRR17181025",
"SRR17181026",
"SRR17181027",
"SRR17181028",
"SRR17181029",
"SRR17181030",
"SRR17181031",
"SRR17181032",
"SRR17391683",
"SRR17391684",
"SRR17391685",
"SRR17391686",
"SRR17573933",
"SRR17573934",
"SRR17573935",
"SRR17573936",
"SRR17573937",
"SRR17573938",
"SRR17573939",
"SRR17573940",
"SRR17573941",
"SRR17573942",
"SRR17573943",
"SRR17573944",
"SRR17573945",
"SRR17573946",
"SRR17573947",
"SRR17573948",
"SRR17573949",
"SRR17573950",
"SRR17573951",
"SRR17573952",
"SRR17573953",
"SRR17573954",
"SRR17573955",
"SRR17573956",
"SRR17573957",
"SRR17573958",
"SRR17573959",
"SRR17573960",
"SRR17573961",
"SRR17573962",
"SRR17573963",
"SRR17573964",
"SRR17573965",
"SRR17573966",
"SRR17573967",
"SRR17573968",
"SRR17573969",
"SRR17573970",
"SRR17573971",
"SRR17573972",
"SRR17573973",
"SRR17573974",
"SRR17573975",
"SRR17573976",
"SRR17573977",
"SRR17573978",
"SRR17573979",
"SRR17573980",
"SRR17573981",
"SRR17573982",
"SRR17573983",
"SRR17573984",
"SRR17573985",
"SRR17573986",
"SRR17573987",
"SRR17573988",
"SRR17573989",
"SRR17573990",
"SRR17573991",
"SRR17573992",
"SRR17573993",
"SRR17573994",
"SRR17573995",
"SRR17573996",
"SRR17573997",
"SRR17573998",
"SRR17573999",
"SRR17574000",
"SRR17574001",
"SRR17574002",
"SRR17574003",
"SRR17574004",
"SRR17574005",
"SRR17574006",
"SRR17574007",
"SRR17574008",
"SRR17574009",
"SRR17574010",
"SRR17574011",
"SRR17574012",
"SRR17574013",
"SRR17574014",
"SRR17574015",
"SRR17574016",
"SRR17574017",
"SRR17574018",
"SRR17574019",
"SRR17574020",
"SRR17574021",
"SRR17574022",
"SRR17574023",
"SRR17574024",
"SRR17574025",
"SRR17574026",
"SRR17574027",
"SRR17574028",
"SRR17574029",
"SRR17574030",
"SRR17574031",
"SRR17574032",
"SRR17574033",
"SRR17574034",
"SRR17574035",
"SRR17574036",
"SRR17574037",
"SRR17574038",
"SRR17574039",
"SRR17574040",
"SRR17574041",
"SRR17574042",
"SRR17574043",
"SRR17574044",
"SRR17574045",
"SRR17574046",
"SRR17574047",
"SRR17574048",
"SRR17574049",
"SRR17574050",
"SRR17574051",
"SRR17574052",
"SRR17574053",
"SRR17574054",
"SRR17574055",
"SRR17574056",
"SRR17574057",
"SRR17574058",
"SRR17574059",
"SRR17574060",
"SRR17574061",
"SRR17574062",
"SRR17574063",
"SRR17574064",
"SRR17574065",
"SRR17574066",
"SRR17574067",
"SRR17574068",
"SRR17574069",
"SRR17574070",
"SRR17574071",
"SRR17574072",
"SRR17574073",
"SRR17574074",
"SRR17574075",
"SRR17574076",
"SRR17574077",
"SRR17574078",
"SRR17574079",
"SRR17574080",
"SRR17574081",
"SRR17574082",
"SRR17574083",
"SRR17574084",
"SRR17574085",
"SRR17574086",
"SRR17574087",
"SRR17574088",
"SRR17574089",
"SRR17574090",
"SRR17574091",
"SRR17574092",
"SRR17574093",
"SRR17574094",
"SRR17574095",
"SRR17574096",
"SRR17574097",
"SRR17574098",
"SRR17574099",
"SRR17574100",
"SRR17574101",
"SRR17574102",
"SRR17574103",
"SRR17574104",
"SRR17574105",
"SRR17574106",
"SRR17574107",
"SRR17574108",
"SRR17574109",
"SRR17574110",
"SRR17574111",
"SRR17574112",
"SRR17574113",
"SRR17574114",
"SRR17574115",
"SRR17574116",
"SRR17574117",
"SRR17574118",
"SRR17574119",
"SRR17574120",
"SRR17574121",
"SRR17574122",
"SRR17574123",
"SRR17574124",
"SRR17574125",
"SRR17574126",
"SRR17574127",
"SRR17574128",
"SRR17574129",
"SRR17660954",
"SRR17660955",
"SRR17967721",
"SRR17970479",
"SRR17970480",
"SRR17970481",
"SRR17970482",
"SRR17970483",
"SRR17970484",
"SRR17970485",
"SRR17970486",
"SRR17971426",
"SRR17971427",
"SRR17971428",
"SRR189663",
"SRR189664",
"SRR189665",
"SRR20117960",
"SRR20117961",
"SRR20117962",
"SRR20117963",
"SRR20117964",
"SRR20117965",
"SRR20117966",
"SRR20117967",
"SRR20117968",
"SRR20117969",
"SRR20117970",
"SRR20117971",
"SRR20117972",
"SRR20117973",
"SRR20117974",
"SRR20117975",
"SRR20117976",
"SRR20117977",
"SRR20117978",
"SRR20117979",
"SRR20117980",
"SRR20117981",
"SRR20117982",
"SRR20117983",
"SRR2954803",
"SRR332251",
"SRR332274",
"SRR332407",
"SRR334209",
"SRR343134",
"SRR343135",
"SRR343136",
"SRR343137",
"SRR343138",
"SRR343139",
"SRR343140",
"SRR343141",
"SRR343142",
"SRR343143",
"SRR343145",
"SRR343146",
"SRR343149",
"SRR343150",
"SRR343152",
"SRR4002443",
"SRR4002444",
"SRR4142426",
"SRR4142434",
"SRR4142436",
"SRR4142462",
"SRR4142832",
"SRR4142833",
"SRR4142834",
"SRR4142839",
"SRR4142844",
"SRR4151122",
"SRR4235915",
"SRR4235916",
"SRR4235917",
"SRR4235920",
"SRR4235921",
"SRR4235922",
"SRR4235923",
"SRR4235924",
"SRR4235925",
"SRR4235926",
"SRR4235927",
"SRR4235929",
"SRR4235930",
"SRR4235931",
"SRR4235932",
"SRR4235933",
"SRR5676587",
"SRR5676590",
"SRR5676591",
"SRR5676592",
"SRR5676593",
"SRR617720",
"SRR617721",
"SRR617722",
"SRR617723",
"SRR617724",
"SRR617725",
"SRR617726",
"SRR617727",
"SRR617728",
"SRR617729",
"SRR617730",
"SRR617731",
"SRR617732",
"SRR617733",
"SRR617734",
"SRR617735",
"SRR617736",
"SRR617737",
"SRR617738",
"SRR617739",
"SRR617740",
"SRR617741",
"SRR617742",
"SRR617743",
"SRR617744",
"SRR617745",
"SRR6277228",
"SRR6277229",
"SRR6434917",
"SRR6434918",
"SRR7418922",
"SRR7418923",
"SRR7418924",
"SRR7418925",
"SRR7418926",
"SRR7418927",
"SRR7418928",
"SRR7418929",
"SRR7418930",
"SRR7418931",
"SRR7418932",
"SRR7418934",
"SRR7418935",
"SRR7418936",
"SRR7418937",
"SRR7418938",
"SRR7418939",
"SRR7418940",
"SRR7418941",
"SRR7418942",
"SRR7418943",
"SRR7418944",
"SRR7418945",
"SRR7418946",
"SRR7418947",
"SRR7418948",
"SRR7418949",
"SRR8759697",
"SRR8759699",
"SRR8759700",
"SRR8759701",
"SRR8759702",
"SRR8759703",
"SRR8759704",
"SRR8759705",
"SRR8759706",
"SRR8759707",
"SRR8759708",
"SRR8759709",
"SRR9067511",
"SRR9265315",
"SRR9265316",
"SRR9265317",
"SRR9265318");
$len=length($hash{$arr[0]});
my @a=();
my @b=();
open O,"veronica/Documents/School/McMaster/PhDBioinformtics/distance_matrix_no_root.txt" or die "Can't open file '$file' for read. $!";
print O "\t";
for $k(0..1859 ){print O $arr[$k],"\t";}
print O "\n";
for $i(0..1859 ){
    for $j($i..1859 ){
        @a=split(undef, $hash{$arr[$i]});
        @b=split(undef, $hash{$arr[$j]});
        #print O $arr[$j],"\t";
        $count=0;
        for $p(0..$len-1){
            if($a[$p] ne $b[$p]){
                $count=$count+1;
            }
        }push @brr, $count;
    }my $str = join "\t", @brr;
    my @brrp=('x') x $i;
    $s = join "\t", @brrp;
    @brr = ();
    if($i==0){print O $arr[$i],"\t",$str,"\n";}else{
        print O $arr[$i],"\t",$s,"\t",$str,"\n";}
    my $str='';
    @brrp = ();
}
close(O);
