//
// Created by layne on 18-9-29.
//

#ifndef ADP_FILTER_COFF_H
#define ADP_FILTER_COFF_H

static constexpr int adp_filter_coeff_len = 1536;
static float adp_filter_coeff[adp_filter_coeff_len] = {
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    -0.1912227, 0.1850990,  -0.1741332, 0.0426635,  -0.2867716, 0.1073588,
    -0.1847986, 0.0758629,  -0.1033363, 0.0989619,  0.0860790,  0.0710217,
    -0.0302731, -0.0417254, -0.0290115, 0.0521865,  0.0461163,  0.0222350,
    -0.1945265, -0.0980438, -0.0664367, 0.1085399,  -0.1671796, 0.1350162,
    0.0321817,  0.0141061,  0.0503558,  -0.0111665, -0.0927427, 0.2204267,
    -0.0614321, -0.0176980, 0.0127379,  0.0205047,  -0.0174790, 0.0531357,
    0.0005383,  0.0063760,  0.0617795,  -0.0085592, 0.0997543,  0.0451434,
    -0.0008809, 0.0142534,  -0.0053794, -0.0192361, 0.0228799,  0.0071443,
    0.0225062,  -0.0000347, 0.0402807,  0.0047772,  0.0637307,  -0.0763051,
    -0.0567031, 0.0340841,  -0.0003984, 0.0073318,  0.0718079,  -0.1146350,
    0.0016296,  -0.0226978, 0.0013946,  -0.0001508, 0.0198202,  -0.0892877,
    -0.0395722, 0.0061120,  0.0186546,  -0.0438137, -0.0615001, -0.0701819,
    -0.0331148, 0.0185835,  -0.0267962, -0.0379473, -0.0286128, -0.0799521,
    -0.0495693, 0.0442766,  -0.0016581, -0.0561241, -0.1613207, 0.0332342,
    -0.0639983, 0.0158932,  -0.1044257, 0.0313127,  -0.1249584, 0.1643697,
    0.0019491,  0.0263484,  0.0081167,  -0.0219748, -0.0439887, 0.2540214,
    0.0218324,  -0.0126702, 0.0089399,  0.0310132,  0.0353820,  0.2717330,
    -0.0071274, -0.0208603, -0.0183486, -0.0714515, 0.0857269,  0.2273994,
    0.0303648,  -0.0299926, 0.0276200,  0.0394801,  0.2037455,  0.1289178,
    0.0195112,  0.0090709,  -0.0344237, -0.0005941, 0.2468869,  0.0159960,
    0.0080312,  -0.0053673, 0.0299719,  -0.0194422, 0.2058839,  -0.0704684,
    0.0187604,  -0.0261185, -0.0296956, 0.0669901,  0.1353050,  -0.1786724,
    -0.0139151, -0.0298198, -0.0117111, -0.0341516, 0.0203791,  -0.2594086,
    0.0076220,  -0.0309799, 0.0189111,  0.0486371,  -0.0633798, -0.2574181,
    0.0150142,  0.0017636,  -0.0206200, -0.0404609, -0.1322614, -0.1709141,
    -0.0261150, -0.0012544, 0.0180059,  0.0003683,  -0.1879186, -0.1143154,
    -0.0307042, -0.0011097, -0.0237162, 0.0031318,  -0.1998704, -0.0296004,
    -0.0114054, 0.0168930,  0.0368511,  -0.0149782, -0.1315464, 0.0491973,
    0.0070868,  0.0363451,  -0.0401009, -0.0022012, -0.0356349, 0.1137286,
    -0.0092395, 0.0270653,  -0.0298016, -0.0317993, -0.0170869, 0.1208340,
    0.0083643,  0.0142911,  0.0068286,  0.0212904,  0.0566630,  0.0899578,
    0.0096766,  -0.0027015, -0.0147317, 0.0012714,  0.0790973,  0.0296583,
    -0.0000523, 0.0218956,  0.0200856,  -0.0039083, 0.1172198,  -0.0056004,
    0.0091093,  -0.0091683, -0.0072361, 0.0167596,  0.0472444,  -0.0572819,
    0.0105002,  -0.0221307, -0.0258911, -0.0113840, -0.0137951, -0.0297111,
    0.0001119,  0.0014754,  0.0089184,  -0.0057895, -0.0267205, -0.0293092,
    0.0214262,  -0.0113195, -0.0176572, 0.0112086,  -0.0522684, -0.0120804,
    0.0119174,  0.0120257,  -0.0127692, -0.0094866, -0.0201765, 0.0197279,
    -0.0058725, 0.0122144,  0.0116927,  0.0062827,  -0.0232587, 0.0161710,
    0.0281816,  -0.0176916, 0.0425937,  0.0051056,  0.0399789,  0.0017286,
    0.0171906,  -0.0319390, 0.0047833,  -0.0048737, 0.0933498,  -0.0569157,
    -0.0168328, 0.0123031,  0.0016354,  -0.0342300, 0.0594845,  -0.0736607,
    -0.0073404, -0.0283699, 0.0028477,  0.0085786,  0.0301351,  -0.1249481,
    0.0079410,  -0.0144183, -0.0008594, 0.0231845,  -0.0050326, -0.1167226,
    -0.0094527, -0.0179354, -0.0168556, -0.0151987, -0.0842205, -0.1326519,
    -0.0115417, -0.0018996, 0.0264984,  0.0101858,  -0.1249921, -0.0754608,
    -0.0181227, -0.0165601, -0.0449608, -0.0334868, -0.0845052, 0.0459656,
    0.0066570,  0.0028092,  0.0046890,  -0.0097163, -0.1229260, 0.1214025,
    0.0084613,  0.0382114,  -0.0318514, 0.0473405,  -0.0505081, 0.1776914,
    0.0110357,  -0.0304949, 0.0206951,  -0.0507101, 0.1209544,  0.1943704,
    0.0226930,  -0.0159607, 0.0393887,  0.0393947,  0.1983817,  0.1404105,
    0.0243941,  0.0204678,  -0.0755725, -0.0060075, 0.2883647,  0.0923487,
    0.0331095,  0.0216956,  0.0724098,  -0.0151255, 0.3173168,  -0.0751199,
    0.0011543,  0.0066405,  -0.0415383, 0.0385291,  0.2254468,  -0.2452478,
    -0.0020078, -0.0040895, 0.0580993,  -0.0582174, 0.1425291,  -0.3895248,
    -0.0051734, 0.0152085,  -0.0598897, 0.0640863,  -0.0210195, -0.4491813,
    0.0212130,  -0.0326879, 0.0008210,  -0.0681873, -0.2499828, -0.4412543,
    -0.0191645, -0.0269455, 0.0674640,  0.0649391,  -0.4879939, -0.3866123,
    -0.0297045, -0.0404516, -0.0764318, -0.0273745, -0.6098565, -0.1559792,
    -0.0335690, -0.0061413, 0.0652716,  -0.0152319, -0.7067572, 0.0653185,
    -0.0159850, -0.0170798, -0.0105365, -0.0010161, -0.6650386, 0.2343277,
    -0.0122895, 0.0323860,  0.0034390,  -0.0864187, -0.5183882, 0.3255857,
    0.0303964,  0.0277773,  -0.0238697, 0.0833965,  -0.3002528, 0.3921800,
    0.0103692,  0.0563071,  -0.0050652, -0.0583882, -0.1194812, 0.3983725,
    0.0143661,  0.0207455,  0.0228823,  0.0841177,  0.0824434,  0.2972976,
    0.0225534,  0.0190451,  -0.0393296, -0.0188309, 0.1573133,  0.2043224,
    0.0095915,  0.0286009,  0.0455019,  -0.0301556, 0.1797914,  0.0075724,
    -0.0013811, -0.0028569, -0.0223396, 0.0237965,  0.1547433,  -0.0533382,
    0.0134802,  0.0092920,  0.0162509,  -0.0167398, 0.0896880,  -0.1411414,
    0.0011939,  -0.0134788, 0.0126892,  0.0511353,  -0.0261814, -0.1479980,
    -0.0043073, -0.0153561, -0.0071907, -0.0304735, -0.0917936, -0.0846464,
    -0.0106143, -0.0011195, 0.0250284,  0.0221013,  -0.1212187, -0.0101443,
    -0.0119129, 0.0087039,  -0.0251466, 0.0041552,  -0.1347371, 0.0722109,
    -0.0112266, 0.0053454,  0.0062918,  -0.0170164, -0.0977245, 0.1089381,
    0.0037313,  0.0123563,  0.0002499,  0.0216194,  -0.0032074, 0.1511802,
    0.0011208,  0.0115708,  -0.0161608, -0.0184872, 0.0548925,  0.1365549,
    0.0029103,  0.0078232,  0.0199025,  0.0275247,  0.1065044,  0.0595130,
    -0.0056795, -0.0005415, -0.0283328, -0.0014079, 0.1346409,  -0.0070070,
    0.0005084,  -0.0078988, 0.0362030,  0.0044571,  0.1180470,  -0.0905686,
    0.0250240,  -0.0015898, -0.0067040, 0.0237563,  0.0510270,  -0.1281550,
    -0.0053828, 0.0044918,  -0.0017161, -0.0237474, -0.0091922, -0.1398301,
    -0.0152658, -0.0030085, 0.0130558,  0.0186122,  -0.0757106, -0.1124480,
    -0.0082354, 0.0015344,  -0.0298202, -0.0087455, -0.1253676, -0.0439571,
    -0.0041849, 0.0054733,  0.0276890,  -0.0037804, -0.1251611, 0.0271266,
    -0.0068049, 0.0083301,  -0.0133083, 0.0226986,  -0.1092472, 0.0813657,
    0.0073594,  0.0106949,  0.0034994,  -0.0149816, -0.0526728, 0.1115131,
    0.0080494,  0.0144086,  0.0039802,  0.0216952,  0.0173442,  0.1063743,
    0.0036843,  -0.0018122, -0.0108216, -0.0212396, 0.0690316,  0.0652805,
    0.0060684,  -0.0010961, 0.0182200,  0.0087184,  0.0947704,  0.0190382,
    -0.0015476, -0.0179560, -0.0207924, 0.0201129,  0.1092614,  -0.0223490,
    -0.0074268, -0.0049806, 0.0095004,  -0.0179828, 0.0764254,  -0.0788782,
    -0.0088571, -0.0036243, 0.0048207,  0.0219131,  0.0187768,  -0.0869143,
    -0.0005726, -0.0053678, -0.0020927, -0.0157797, -0.0313469, -0.0909789,
    -0.0082654, -0.0022259, 0.0210155,  0.0107800,  -0.0657631, -0.0569130,
    -0.0051553, 0.0005189,  -0.0118931, -0.0052730, -0.0730304, -0.0135517,
    -0.0031049, -0.0042033, 0.0066737,  -0.0047335, -0.0621953, 0.0257897,
    0.0006342,  0.0041807,  -0.0062111, 0.0120581,  -0.0216699, 0.0673214,
    -0.0013644, 0.0036253,  -0.0012042, -0.0081968, 0.0199624,  0.0624616,
    0.0054624,  0.0018992,  0.0066751,  0.0084560,  0.0490699,  0.0412383,
    0.0004272,  -0.0051843, -0.0119557, 0.0013947,  0.0624459,  0.0132529,
    0.0015621,  -0.0126005, 0.0150151,  -0.0071002, 0.0583337,  -0.0188580,
    0.0021251,  -0.0040499, -0.0022039, 0.0078918,  0.0355540,  -0.0496256,
    -0.0036419, -0.0017119, -0.0021250, -0.0062566, -0.0011109, -0.0531303,
    -0.0012022, -0.0016638, 0.0066830,  0.0051417,  -0.0172834, -0.0431736,
    -0.0029410, -0.0017266, -0.0076120, -0.0063793, -0.0336428, -0.0230979,
    -0.0041878, 0.0036256,  0.0073266,  -0.0029472, -0.0337047, -0.0031642,
    -0.0026587, -0.0039698, -0.0091414, -0.0013359, -0.0279975, 0.0143442,
    -0.0474144, 0.0060696,  -0.0296154, 0.0006136,  -0.0068833, 0.0705231,
    0.0320243,  0.0257965,  0.0065879,  -0.0264677, -0.0702915, 0.1442964,
    -0.0095345, -0.0312792, 0.0416841,  -0.0152720, -0.0551221, -0.1396452,
    0.0019089,  -0.0402145, -0.0698273, 0.0745254,  -0.0300285, -0.0378171,
    -0.1344261, -0.0656269, -0.0077642, -0.0383189, -0.0852450, -0.0298054,
    -0.2571020, 0.0365747,  0.1669619,  -0.1468360, 0.0975295,  -0.3010011,
    -0.0727764, -0.2137894, 0.1038804,  0.1982930,  0.0333594,  -0.1747120,
    0.0183768,  -0.2087581, -0.0188623, 0.1927662,  -0.0769083, -0.0987765,
    0.0606294,  0.0061882,  0.0146982,  0.0362168,  0.2690670,  0.0030640,
    0.0085473,  0.1074151,  -0.0640013, 0.0175161,  -0.0582388, -0.0405226,
    -0.0852050, -0.0048755, 0.0185345,  0.0956255,  0.0226524,  -0.0383554,
    0.0718989,  -0.0369001, -0.0474142, 0.0105181,  -0.1016836, -0.0218953,
    -0.0274312, -0.0360292, 0.0349053,  -0.0304307, -0.0445732, 0.0751726,
    0.0104603,  0.0235997,  0.0097164,  -0.0056312, 0.0381956,  -0.0163091,
    -0.0003253, -0.0134541, 0.0093706,  0.0324109,  -0.0546963, 0.0306905,
    0.0514693,  -0.0163546, -0.0087376, -0.0424844, 0.0118898,  0.0164375,
    0.0114150,  0.0072139,  -0.0070355, 0.0034848,  0.0091579,  0.0032901,
    -0.0027745, -0.0005701, 0.0033587,  0.0016355,  -0.0036181, -0.0105640,
    -0.0015214, 0.0006795,  -0.0031790, -0.0020504, 0.0067836,  -0.0017227,
    0.0041087,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
    0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000};
#endif  // ADP_FILTER_COFF_H
