{-# OPTIONS_GHC -Wall #-}
module Anno where

import Data.List

-- What does this header mean?
-- >Feature ref|NC_012920.1|

to_tab :: String -> [Anno] -> [String]
to_tab nm ann = (">Feature " ++ nm) : (map (intercalate "\t") $ concatMap to_tab1 ann)
  where
    to_tab1 :: Anno -> [[String]]
    to_tab1 (Gene s e n syns what ps ns) =
        [ show s, show e, label what n ] :
        ( if has_gene what then [[ "", "", "", "gene", n ]] else [] ) ++
        [ [ "", "", "", "gene_syn", sy ] | sy <- syns ] ++
        [ [ show s, show e, w ] | w <- describe what ] ++
        [ [ "", "", "", "product", p ] | p <- [ps], not (null p) ] ++
        [ [ "", "", "", "note", n ] | n <- ns ] ++
        more what

    describe CDS = ["CDS"]
    describe CDS' = ["CDS"]
    describe TRNA = ["tRNA"]
    describe RRNA = ["rRNA"]
    describe _ = []

    has_gene (STS _) = False
    has_gene Other = False
    has_gene _ = True

    label (STS _) _ = "STS"
    label Other n = n
    label _ _ = "gene"

    more CDS' = [ [ "", "", "", "note", "TAA stop codon is completed by the addition of 3' A residues to the mRNA" ] ]
    more (STS sn) = [ [ "", "", "", "standard_name", sn ] ]
    more _ = []

data Anno
    = Gene { start :: Int
           , end :: Int
           , name :: String
           , syn :: [String]
           , what :: What
           , product :: String
           , note :: [String] }
    deriving Show

data What = CDS | CDS' | TRNA | RRNA | Other | STS String deriving Show

rCRS_anno :: [Anno]
rCRS_anno =
    [ Gene 576	1	"D-loop" [] Other "" []
    , Gene 16569	16024 "" [] Other "" []
    , Gene 577	647	    "TRNF" [] TRNA "tRNA-Phe" []
    , Gene 648	1601	"RNR1" ["MTRNR1"] RRNA "s-rRNA" ["12S rRNA; 12S ribosomal RNA"]
    , Gene 1602	1670	"TRNV" [] TRNA "tRNA-Val" []
    , Gene 1671	3229	"RNR2" [] RRNA "l-rRNA" ["16S ribosomal RNA; 16S rRNA"]
    , Gene 3230	3304   "TRNL1" ["MTTL1"] TRNA "tRNA-Leu" []
    , Gene 3307	4262	"ND1"  [] CDS' "NADH dehydrogenase subunit 1" []

    , Gene 4263	4331    "TRNI" [] TRNA "tRNA-Ile" []
    , Gene 4400	4329	"TRNQ" [] TRNA "tRNA-Gln" []
    , Gene 4402	4469	"TRNM" [] TRNA "tRNA-Met" []
    , Gene 4470	5511	"ND2"  [] CDS' "NADH dehydrogenase subunit 2" []
    , Gene 5512	5579	"TRNW" [] TRNA "tRNA-Trp" []
    , Gene 5655	5587	"TRNA" [] TRNA "tRNA-Ala" []
    , Gene 5729	5657	"TRNN" [] TRNA "tRNA-Asn" []
    , Gene 5826	5761	"TRNC" [] TRNA "tRNA-Cys" []
    , Gene 5891	5826	"TRNY" [] TRNA "tRNA-Tyr" []
    , Gene 5904	7445	"COX1" ["COI"] CDS "cytochrome c oxidase subunit I" ["cytochrome c oxidase I"]
    , Gene 7514	7446   "TRNS1" [] TRNA "tRNA-Ser" []
    , Gene 7518	7585	"TRND" [] TRNA "tRNA-Asp" []
    , Gene 7586	8269	"COX2" [] CDS "cytochrome c oxidase subunit II" ["cytochrome c oxidase II"]
    , Gene 8295	8364	"TRNK" [] TRNA "tRNA-Lys" []
    , Gene 8366	8572	"ATP8" [] CDS  "ATP synthase F0 subunit 8" ["ATP synthase 8; ATPase subunit 8"]
    , Gene 8527	9207	"ATP6" [] CDS  "ATP synthase F0 subunit 6" ["ATP synthase 6; ATPase subunit 6"]
    , Gene 9207	9990	"COX3" [] CDS' "cytochrome c oxidase subunit III" []
    , Gene 9342	9416	""     [] (STS "PMC55343P8") "" []
    , Gene 9991	    10058	"TRNG" [] TRNA "tRNA-Gly" []
    , Gene 10059	10404	"ND3"  [] CDS' "NADH dehydrogenase subunit 3" []
    , Gene 10405	10469	"TRNR" [] TRNA "tRNA-Arg" []
    , Gene 10470	10766	"ND4L" [] CDS "NADH dehydrogenase subunit 4L" []
    , Gene 10760	12137	"ND4"  [] CDS' "NADH dehydrogenase subunit 4" []
    , Gene 12138	12206	"TRNH" [] TRNA "tRNA-His" []
    , Gene 12207	12265  "TRNS2" [] TRNA "tRNA-Ser" []
    , Gene 12266	12336  "TRNL2" [] TRNA "tRNA-Leu" []
    , Gene 12337	14148	"ND5"  [] CDS  "NADH dehydrogenase subunit 5" []
    , Gene 14673	14149	"ND6"  [] CDS  "NADH dehydrogenase subunit 6" []
    , Gene 14742	14674	"TRNE" [] TRNA "tRNA-Glu" []
    , Gene 14747	15887	"CYTB" [] CDS' "cytochrome b" []
    , Gene 15888	15953	"TRNT" [] TRNA "tRNA-Thr" []
    , Gene 16023	15956	"TRNP" [] TRNA "tRNA-Pro" [] ]

aas :: [(String, String)]
aas = [
     (,) "ND1"
                     "MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYG\
                     \LLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGL\
                     \LFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYEVTLAIILLSTLLMSGSF\
                     \NLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAG\
                     \PFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTA\
                     \YPRFRYDQLMHLLWKNFLPLTLALLMWYVSMPITISSIPPQT",
     (,) "ND2"
                     "MNPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVL\
                     \TKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAM\
                     \AMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPSLNVSLLLTLSIL\
                     \SIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLL\
                     \LNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSL\
                     \IIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKMKWQFEHTKPTPFLPTLIALTT\
                     \LLLPISPFMLMIL",
     (,) "COX1"
                     "MFADRWLFSTNHKDIGTLYLLFGAWAGVLGTALSLLIRAELGQP\
                     \GNLLGNDHIYNVIVTAHAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSF\
                     \WLLPPSLLLLLASAMVEAGAGTGWTVYPPLAGNYSHPGASVDLTIFSLHLAGVSSILG\
                     \AINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT\
                     \FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIVTYYSGKKEPFGYMGMVWA\
                     \MMSIGFLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGSNMKW\
                     \SAAVLWALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGG\
                     \FIHWFPLFSGYTLDQTYAKIHFTIMFIGVNLTFFPQHFLGLSGMPRRYSDYPDAYTTW\
                     \NILSSVGSFISLTAVMLMIFMIWEAFASKRKVLMVEEPSMNLEWLYGCPPPYHTFEEP\
                     \VYMKS",
     (,) "COX2"
                     "MAHAAQVGLQDATSPIMEELITFHDHALMIIFLICFLVLYALFL\
                     \TLTTKLTNTNISDAQEMETVWTILPAIILVLIALPSLRILYMTDEVNDPSLTIKSIGH\
                     \QWYWTYEYTDYGGLIFNSYMLPPLFLEPGDLRLLDVDNRVVLPIEAPIRMMITSQDVL\
                     \HSWAVPTLGLKTDAIPGRLNQTTFTATRPGVYYGQCSEICGANHSFMPIVLELIPLKI\
                     \FEMGPVFTL",
     (,) "ATP8"
                     "MPQLNTTVWPTMITPMLLTLFLITQLKMLNTNYHLPPSPKPMKM\
                     \KNYNKPWEPKWTKICSLHSLPPQS",
     (,) "ATP6"
                     "MNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLIT\
                     \TQQWLIKLTSKQMMTMHNTKGRTWSLMLVSLIIFIATTNLLGLLPHSFTPTTQLSMNL\
                     \AMAIPLWAGTVIMGFRSKIKNALAHFLPQGTPTPLIPMLVIIETISLLIQPMALAVRL\
                     \TANITAGHLLMHLIGSATLAMSTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLV\
                     \SLYLHDNT",
     (,) "COX3"
                     "MTHQSHAYHMVKPSPWPLTGALSALLMTSGLAMWFHFHSMTLLM\
                     \LGLLTNTLTMYQWWRDVTRESTYQGHHTPPVQKGLRYGMILFITSEVFFFAGFFWAFY\
                     \HSSLAPTPQLGGHWPPTGITPLNPLEVPLLNTSVLLASGVSITWAHHSLMENNRNQMI\
                     \QALLITILLGLYFTLLQASEYFESPFTISDGIYGSTFFVATGFHGLHVIIGSTFLTIC\
                     \FIRQLMFHFTSKHHFGFEAAAWYWHFVDVVWLFLYVSIYWWGS",
     (,) "ND3"
                     "MNFALILMINTLLALLLMIITFWLPQLNGYMEKSTPYECGFDPM\
                     \SPARVPFSMKFFLVAITFLLFDLEIALLLPLPWALQTTNLPLMVMSSLLLIIILALSL\
                     \AYEWLQKGLDWTE",
     (,) "ND4L"
                     "MPLIYMNIMLAFTISLLGMLVYRSHLMSSLLCLEGMMLSLFIMA\
                     \TLMTLNTHSLLANIVPIAMLVFAACEAAVGLALLVSISNTYGLDYVHNLNLLQC",
     (,) "ND4"
                     "MLKLIVPTIMLLPLTWLSKKHMIWINTTTHSLIISIIPLLFFNQ\
                     \INNNLFSCSPTFSSDPLTTPLLMLTTWLLPLTIMASQRHLSSEPLSRKKLYLSMLISL\
                     \QISLIMTFTATELIMFYIFFETTLIPTLAIITRWGNQPERLNAGTYFLFYTLVGSLPL\
                     \LIALIYTHNTLGSLNILLLTLTAQELSNSWANNLMWLAYTMAFMVKMPLYGLHLWLPK\
                     \AHVEAPIAGSMVLAAVLLKLGGYGMMRLTLILNPLTKHMAYPFLVLSLWGMIMTSSIC\
                     \LRQTDLKSLIAYSSISHMALVVTAILIQTPWSFTGAVILMIAHGLTSSLLFCLANSNY\
                     \ERTHSRIMILSQGLQTLLPLMAFWWLLASLANLALPPTINLLGELSVLVTTFSWSNIT\
                     \LLLTGLNMLVTALYSLYMFTTTQWGSLTHHINNMKPSFTRENTLMFMHLSPILLLSLN\
                     \PDIITGFSS",
     (,) "ND5"
                     "MTMHTTMTTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTF\
                     \IISLFPTTMFMCLDQEVIISNWHWATTQTTQLSLSFKLDYFSMMFIPVALFVTWSIME\
                     \FSLWYMNSDPNINQFFKYLLIFLITMLILVTANNLFQLFIGWEGVGIMSFLLISWWYA\
                     \RADANTAAIQAILYNRIGDIGFILALAWFILHSNSWDPQQMALLNANPSLTPLLGLLL\
                     \AAAGKSAQLGLHPWLPSAMEGPTPVSALLHSSTMVVAGIFLLIRFHPLAENSPLIQTL\
                     \TLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLMMVTIGINQPHLAFLHICTHAF\
                     \FKAMLFMCSGSIIHNLNNEQDIRKMGGLLKTMPLTSTSLTIGSLALAGMPFLTGFYSK\
                     \DHIIETANMSYTNAWALSITLIATSLTSAYSTRMILLTLTGQPRFPTLTNINENNPTL\
                     \LNPIKRLAAGSLFAGFLITNNISPASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLT\
                     \NKLKMKSPLCTFYFSNMLGFYPSITHRTIPYLGLLTSQNLPLLLLDLTWLEKLLPKTI\
                     \SQHQISTSIITSTQKGMIKLYFLSFFFPLILTLLLIT",
     (,) "ND6"
                     "MMYALFLLSVGLVMGFVGFSSKPSPIYGGLVLIVSGVVGCVIIL\
                     \NFGGGYMGLMVFLIYLGGMMVVFGYTTAMAIEEYPEAWGSGVEVLVSVLVGLAMEVGL\
                     \VLWVKEYDGVVVVVNFNSVGSWMIYEGEGSGLIREDPIGAGALYDYGRWLVVVTGWTL\
                     \FVGVYIVIEIARGN",
     (,) "CYTB"
                     "MTPMRKTNPLMKLINHSFIDLPTPSNISAWWNFGSLLGACLILQ\
                     \ITTGLFLAMHYSPDASTAFSSIAHITRDVNYGWIIRYLHANGASMFFICLFLHIGRGL\
                     \YYGSFLYSETWNIGIILLLATMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTDL\
                     \VQWIWGGYSVDSPTLTRFFTFHFILPFIIAALATLHLLFLHETGSNNPLGITSHSDKI\
                     \TFHPYYTIKDALGLLLFLLSLMTLTLFSPDLLGDPDNYTLANPLNTPPHIKPEWYFLF\
                     \AYTILRSVPNKLGGVLALLLSILILAMIPILHMSKQQSMMFRPLSQSLYWLLAADLLI\
                     \LTWIGGQPVSYPFTIIGQVASVLYFTTILILMPTISLIENKMLKWA" ]
