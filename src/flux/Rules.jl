function calculate_rules_vector(data_dictionary::Dict{String,Any}, value_array::Array{Float64,1})

	# Initialize -
	v = zeros(381)

	function estimate_default(data_dictionary)
		return data_dictionary["default_vmax_value"]
	end
	default() = estimate_default(data_dictionary)

	function calculate_fn_value(x, value_array)
		return value_array[x]
	end
	fn(x) = calculate_fn_value(x, value_array)

	v[1] = default()
	v[2] = default()
	v[3] = default()
	v[4] = default()
	v[5] = default()
	v[6] = default()
	v[7] = default()
	v[8] = default()
	v[9] = default()
	v[10] = default()
	v[11] = default()
	v[12] = default()
	v[13] = default()
	v[14] = default()
	v[15] = default()
	v[16] = default()
	v[17] = default()
	v[18] = default()
	v[19] = default()
	v[20] = default()
	v[21] = default()
	v[22] = default()
	v[23] = default()
	v[24] = default()
	v[25] = default()
	v[26] = default()
	v[27] = default()
	v[28] = default()
	v[29] = default()
	v[30] = default()
	v[31] = default()
	v[32] = default()
	v[33] = default()
	v[34] = fn(1922)
	v[35] = max(min(min(min(fn(346), fn(1764)), fn(347)), fn(1047)), min(min(min(fn(346), fn(1764)), fn(347)), fn(1048)))
	v[36] = fn(1922)
	v[37] = fn(814)
	v[38] = fn(552)
	v[39] = fn(511)
	v[40] = default()
	v[41] = default()
	v[42] = default()
	v[43] = default()
	v[44] = fn(1597)
	v[45] = default()
	v[46] = max(fn(1096), fn(1097))
	v[47] = fn(886)
	v[48] = fn(881)
	v[49] = max(min(fn(769), fn(770)), fn(881))
	v[50] = fn(1578)
	v[51] = max(max(max(max(max(max(max(max(max(fn(791), fn(793)), fn(790)), fn(792)), fn(794)), min(fn(791), fn(813))), min(fn(792), fn(813))), min(fn(793), fn(813))), min(fn(794), fn(813))), min(fn(790), fn(813)))
	v[52] = max(fn(976), fn(977))
	v[53] = max(fn(849), fn(828))
	v[54] = fn(870)
	v[55] = max(fn(828), fn(641))
	v[56] = fn(1049)
	v[57] = fn(1710)
	v[58] = max(max(max(fn(411), fn(608)), fn(607)), fn(177))
	v[59] = max(fn(607), fn(608))
	v[60] = default()
	v[61] = max(fn(286), fn(287))
	v[62] = max(fn(286), fn(287))
	v[63] = max(max(fn(298), fn(179)), fn(178))
	v[64] = max(max(max(max(max(max(fn(1377), fn(1334)), fn(1376)), fn(245)), fn(48)), fn(49)), fn(1375))
	v[65] = fn(987)
	v[66] = fn(59)
	v[67] = max(min(min(min(fn(346), fn(1764)), fn(347)), fn(1047)), min(min(min(fn(346), fn(1764)), fn(347)), fn(1048)))
	v[68] = fn(1801)
	v[69] = max(fn(711), fn(1819))
	v[70] = fn(422)
	v[71] = default()
	v[72] = default()
	v[73] = fn(883)
	v[74] = fn(938)
	v[75] = max(fn(948), fn(949))
	v[76] = fn(128)
	v[77] = max(fn(1792), fn(16))
	v[78] = default()
	v[79] = default()
	v[80] = default()
	v[81] = default()
	v[82] = default()
	v[83] = fn(1752)
	v[84] = max(max(fn(944), fn(945)), fn(946))
	v[85] = max(max(max(fn(128), fn(1321)), fn(1290)), fn(1773))
	v[86] = default()
	v[87] = fn(1711)
	v[88] = max(fn(1851), fn(15))
	v[89] = fn(693)
	v[90] = fn(694)
	v[91] = max(max(max(max(fn(1517), fn(1518)), fn(1519)), fn(1524)), fn(1525))
	v[92] = max(max(fn(724), fn(726)), fn(727))
	v[93] = default()
	v[94] = fn(1647)
	v[95] = fn(1073)
	v[96] = fn(824)
	v[97] = fn(1454)
	v[98] = fn(1073)
	v[99] = fn(1711)
	v[100] = max(fn(236), fn(237))
	v[101] = max(fn(45), fn(44))
	v[102] = fn(39)
	v[103] = max(fn(37), fn(1383))
	v[104] = max(fn(1109), fn(1981))
	v[105] = max(max(max(fn(116), fn(117)), fn(120)), fn(121))
	v[106] = fn(1440)
	v[107] = fn(1933)
	v[108] = fn(1588)
	v[109] = max(max(fn(1514), fn(1515)), fn(1516))
	v[110] = fn(1305)
	v[111] = default()
	v[112] = default()
	v[113] = default()
	v[114] = default()
	v[115] = default()
	v[116] = default()
	v[117] = fn(1710)
	v[118] = default()
	v[119] = fn(246)
	v[120] = max(fn(251), fn(252))
	v[121] = max(fn(1367), fn(1368))
	v[122] = fn(1182)
	v[123] = fn(1456)
	v[124] = fn(1456)
	v[125] = fn(341)
	v[126] = fn(340)
	v[127] = max(fn(1846), fn(1847))
	v[128] = fn(342)
	v[129] = fn(343)
	v[130] = fn(1711)
	v[131] = fn(1429)
	v[132] = fn(458)
	v[133] = fn(966)
	v[134] = default()
	v[135] = fn(382)
	v[136] = fn(65)
	v[137] = max(fn(387), min(fn(769), fn(770)))
	v[138] = max(max(fn(1313), fn(387)), min(fn(769), fn(770)))
	v[139] = max(fn(387), min(fn(769), fn(770)))
	v[140] = max(max(max(max(fn(409), fn(410)), fn(601)), fn(407)), fn(408))
	v[141] = max(min(fn(426), fn(428)), min(fn(426), fn(427)))
	v[142] = fn(431)
	v[143] = max(max(max(max(max(fn(1179), fn(444)), fn(461)), fn(462)), fn(445)), fn(446))
	v[144] = fn(444)
	v[145] = default()
	v[146] = default()
	v[147] = fn(450)
	v[148] = max(max(max(max(fn(467), fn(477)), fn(469)), fn(468)), fn(472))
	v[149] = fn(75)
	v[150] = default()
	v[151] = default()
	v[152] = default()
	v[153] = fn(75)
	v[154] = fn(448)
	v[155] = fn(470)
	v[156] = fn(255)
	v[157] = fn(695)
	v[158] = default()
	v[159] = fn(1419)
	v[160] = fn(555)
	v[161] = max(max(fn(1676), fn(1677)), fn(1675))
	v[162] = max(fn(588), fn(617))
	v[163] = max(fn(596), fn(597))
	v[164] = min(min(min(fn(681), fn(346)), fn(627)), fn(671))
	v[165] = min(min(min(fn(681), fn(346)), fn(627)), fn(671))
	v[166] = min(min(min(fn(681), fn(346)), fn(627)), fn(671))
	v[167] = max(max(max(min(fn(748), fn(762)), min(fn(1916), fn(761))), min(fn(748), fn(761))), min(fn(1916), fn(762)))
	v[168] = max(fn(1477), fn(1478))
	v[169] = fn(744)
	v[170] = fn(616)
	v[171] = max(max(max(max(max(max(max(max(max(max(max(max(fn(269), fn(763)), fn(1371)), fn(1372)), fn(1769)), fn(1599)), fn(262)), fn(1527)), fn(1526)), fn(1528)), fn(1529)), fn(114)), fn(253))
	v[172] = max(max(max(min(fn(748), fn(762)), min(fn(1916), fn(762))), min(fn(748), fn(761))), min(fn(1916), fn(761)))
	v[173] = max(max(max(fn(128), fn(1321)), fn(1290)), fn(1773))
	v[174] = default()
	v[175] = fn(1419)
	v[176] = max(fn(680), fn(679))
	v[177] = max(fn(680), fn(679))
	v[178] = max(max(fn(677), fn(661)), fn(662))
	v[179] = fn(1308)
	v[180] = max(fn(618), fn(619))
	v[181] = max(max(fn(1387), fn(1386)), fn(535))
	v[182] = max(max(max(fn(1851), fn(15)), fn(1731)), fn(1790))
	v[183] = max(max(max(max(fn(1517), fn(1518)), fn(1519)), fn(1524)), fn(1525))
	v[184] = fn(606)
	v[185] = fn(501)
	v[186] = fn(422)
	v[187] = max(max(max(fn(128), fn(1321)), fn(1290)), fn(1773))
	v[188] = default()
	v[189] = fn(1897)
	v[190] = fn(1406)
	v[191] = fn(458)
	v[192] = max(max(max(max(max(max(max(max(max(max(fn(1699), fn(1707)), fn(877)), fn(1700)), fn(1704)), fn(525)), fn(1703)), fn(1708)), fn(1697)), fn(1706)), fn(1698))
	v[193] = default()
	v[194] = fn(1536)
	v[195] = default()
	v[196] = fn(842)
	v[197] = default()
	v[198] = max(max(fn(767), min(fn(769), fn(770))), fn(771))
	v[199] = fn(767)
	v[200] = default()
	v[201] = fn(1987)
	v[202] = max(max(max(max(max(max(max(max(max(max(fn(621), fn(622)), fn(623)), fn(784)), fn(785)), fn(786)), fn(787)), fn(788)), fn(789)), fn(795)), fn(1755))
	v[203] = fn(780)
	v[204] = fn(109)
	v[205] = fn(1917)
	v[206] = fn(801)
	v[207] = max(fn(803), fn(802))
	v[208] = default()
	v[209] = max(fn(800), fn(1293))
	v[210] = max(max(max(max(max(max(fn(1670), fn(1668)), fn(1669)), fn(1671)), fn(1983)), fn(1924)), fn(1925))
	v[211] = default()
	v[212] = default()
	v[213] = max(max(max(max(max(min(min(fn(831), fn(833)), fn(836)), min(min(fn(831), fn(834)), fn(836))), min(min(fn(831), fn(835)), fn(836))), min(min(fn(831), fn(833)), fn(836))), min(min(fn(831), fn(834)), fn(837))), min(min(fn(831), fn(835)), fn(837)))
	v[214] = fn(830)
	v[215] = fn(1428)
	v[216] = fn(128)
	v[217] = fn(987)
	v[218] = max(max(fn(852), fn(853)), fn(854))
	v[219] = max(fn(838), fn(1948))
	v[220] = fn(450)
	v[221] = max(max(max(fn(1403), fn(1573)), fn(2033)), fn(69))
	v[222] = fn(1848)
	v[223] = max(max(max(fn(1578), fn(1950)), fn(516)), fn(1939))
	v[224] = max(max(max(max(max(max(max(fn(890), fn(891)), fn(892)), fn(893)), fn(1955)), fn(305)), fn(1330)), min(fn(890), fn(891)))
	v[225] = fn(1428)
	v[226] = max(max(max(max(max(max(fn(128), fn(1532)), fn(1321)), fn(1773)), fn(1554)), fn(1555)), fn(1825))
	v[227] = fn(303)
	v[228] = fn(901)
	v[229] = fn(899)
	v[230] = fn(1455)
	v[231] = max(fn(128), fn(1596))
	v[232] = max(fn(1792), fn(16))
	v[233] = max(max(fn(732), fn(733)), fn(734))
	v[234] = min(fn(1380), fn(1464))
	v[235] = max(fn(210), fn(922))
	v[236] = fn(923)
	v[237] = fn(924)
	v[238] = max(fn(925), fn(84))
	v[239] = fn(967)
	v[240] = fn(1313)
	v[241] = max(fn(850), fn(851))
	v[242] = fn(1134)
	v[243] = fn(1816)
	v[244] = fn(965)
	v[245] = default()
	v[246] = max(max(fn(955), fn(1746)), fn(707))
	v[247] = max(fn(579), fn(73))
	v[248] = fn(955)
	v[249] = fn(73)
	v[250] = default()
	v[251] = default()
	v[252] = default()
	v[253] = default()
	v[254] = max(max(fn(1538), fn(1540)), fn(1536))
	v[255] = max(max(max(max(fn(1565), fn(1566)), fn(885)), fn(1567)), fn(1564))
	v[256] = fn(382)
	v[257] = max(max(max(max(max(min(fn(1019), fn(1021)), min(fn(1020), fn(1021))), fn(1022)), fn(751)), fn(21)), fn(750))
	v[258] = fn(1023)
	v[259] = max(max(max(max(max(min(fn(1019), fn(1021)), min(fn(1020), fn(1021))), fn(1022)), fn(21)), fn(750)), fn(751))
	v[260] = max(max(max(max(max(min(fn(1019), fn(1021)), min(fn(1020), fn(1021))), fn(1022)), fn(21)), fn(750)), fn(751))
	v[261] = max(max(max(max(max(min(fn(1019), fn(1021)), min(fn(1020), fn(1021))), fn(1022)), fn(21)), fn(750)), fn(751))
	v[262] = max(max(max(max(max(min(fn(1020), fn(1021)), min(fn(1019), fn(1021))), fn(1022)), fn(21)), fn(750)), fn(751))
	v[263] = max(max(max(max(max(min(fn(1019), fn(1021)), min(fn(1020), fn(1021))), fn(1022)), fn(21)), fn(750)), fn(751))
	v[264] = max(max(max(max(max(min(fn(1019), fn(1021)), min(fn(1020), fn(1021))), fn(1022)), fn(21)), fn(750)), fn(751))
	v[265] = max(fn(1570), fn(1571))
	v[266] = max(fn(1570), fn(1571))
	v[267] = max(max(max(max(max(max(max(max(max(max(min(fn(1016), fn(1013)), min(fn(1012), fn(499))), min(fn(1018), fn(1011))), min(fn(1016), fn(1015))), min(fn(1017), fn(1013))), min(fn(1017), fn(1011))), min(fn(1017), fn(1012))), min(fn(1016), fn(1011))), min(fn(499), fn(1011))), min(fn(1015), fn(499))), min(fn(1016), fn(1012)))
	v[268] = default()
	v[269] = default()
	v[270] = default()
	v[271] = default()
	v[272] = fn(1050)
	v[273] = max(min(min(min(fn(1432), fn(1433)), fn(320)), fn(346)), min(min(min(fn(1432), fn(1434)), fn(320)), fn(346)))
	v[274] = max(min(min(min(fn(1432), fn(1433)), fn(320)), fn(346)), min(min(min(fn(1432), fn(1434)), fn(320)), fn(346)))
	v[275] = max(min(min(min(fn(1432), fn(1433)), fn(320)), fn(346)), min(min(min(fn(1432), fn(1434)), fn(320)), fn(346)))
	v[276] = fn(1685)
	v[277] = fn(1041)
	v[278] = max(fn(1792), fn(16))
	v[279] = max(fn(16), fn(1792))
	v[280] = max(max(fn(1554), fn(1555)), fn(1825))
	v[281] = fn(1685)
	v[282] = max(fn(1863), fn(1864))
	v[283] = max(fn(1417), fn(1418))
	v[284] = fn(208)
	v[285] = max(fn(1074), fn(1075))
	v[286] = max(min(min(fn(345), min(fn(346), fn(1764))), min(fn(1174), fn(1175))), min(min(fn(345), min(fn(346), fn(1764))), min(fn(1169), fn(1175))))
	v[287] = default()
	v[288] = max(max(max(max(max(max(max(max(fn(1195), fn(1197)), fn(1196)), min(fn(1195), fn(1197))), min(fn(1196), fn(1197))), fn(1198)), min(fn(1195), fn(1198))), min(fn(1196), fn(1198))), min(fn(1197), fn(1198)))
	v[289] = fn(599)
	v[290] = fn(697)
	v[291] = max(max(fn(1204), fn(1205)), fn(844))
	v[292] = fn(566)
	v[293] = max(max(max(max(fn(1199), fn(1200)), fn(1201)), fn(1607)), fn(1608))
	v[294] = max(fn(1206), fn(1329))
	v[295] = fn(1986)
	v[296] = default()
	v[297] = fn(1066)
	v[298] = max(max(max(max(max(max(fn(128), fn(1321)), fn(1290)), fn(1773)), fn(1554)), fn(1555)), fn(1825))
	v[299] = default()
	v[300] = max(fn(1586), fn(1587))
	v[301] = default()
	v[302] = max(max(fn(1797), fn(1798)), fn(1247))
	v[303] = fn(62)
	v[304] = fn(1301)
	v[305] = max(max(max(fn(1853), fn(1854)), fn(1856)), fn(1855))
	v[306] = min(fn(1079), fn(1080))
	v[307] = default()
	v[308] = default()
	v[309] = max(fn(596), fn(597))
	v[310] = max(fn(596), fn(597))
	v[311] = fn(59)
	v[312] = fn(1187)
	v[313] = fn(422)
	v[314] = max(max(fn(1321), fn(1290)), fn(1773))
	v[315] = default()
	v[316] = max(max(max(fn(1365), fn(1366)), fn(454)), fn(455))
	v[317] = fn(528)
	v[318] = max(fn(759), fn(760))
	v[319] = max(fn(38), fn(1392))
	v[320] = max(fn(1885), fn(1393))
	v[321] = fn(2012)
	v[322] = max(max(max(max(max(fn(1241), fn(1242)), fn(1243)), fn(1244)), fn(1245)), fn(145))
	v[323] = fn(1578)
	v[324] = min(fn(1451), fn(1452))
	v[325] = min(fn(1451), fn(1452))
	v[326] = min(fn(1451), fn(1452))
	v[327] = min(fn(1451), fn(1452))
	v[328] = max(fn(1443), fn(1444))
	v[329] = fn(475)
	v[330] = fn(13)
	v[331] = fn(13)
	v[332] = min(fn(50), fn(1991))
	v[333] = max(max(max(fn(128), fn(1321)), fn(1290)), fn(1773))
	v[334] = fn(584)
	v[335] = fn(1601)
	v[336] = fn(1610)
	v[337] = fn(457)
	v[338] = min(fn(1888), fn(1889))
	v[339] = fn(1633)
	v[340] = max(fn(1076), fn(1077))
	v[341] = max(max(fn(128), fn(1290)), fn(1773))
	v[342] = max(max(fn(1646), fn(1780)), fn(1802))
	v[343] = max(max(fn(1646), fn(1780)), fn(1802))
	v[344] = fn(1665)
	v[345] = max(fn(1649), fn(709))
	v[346] = max(max(max(fn(1661), fn(1662)), fn(1664)), fn(1663))
	v[347] = fn(1637)
	v[348] = fn(128)
	v[349] = fn(693)
	v[350] = fn(128)
	v[351] = fn(1182)
	v[352] = max(max(max(fn(1575), fn(1776)), fn(1540)), fn(1536))
	v[353] = fn(864)
	v[354] = fn(1428)
	v[355] = fn(128)
	v[356] = fn(603)
	v[357] = max(max(fn(1468), fn(1467)), fn(1384))
	v[358] = default()
	v[359] = default()
	v[360] = default()
	v[361] = fn(1202)
	v[362] = default()
	v[363] = default()
	v[364] = default()
	v[365] = default()
	v[366] = max(max(max(max(max(max(max(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(fn(1361), fn(1091)), fn(957)), fn(958)), fn(959)), fn(960)), fn(961)), fn(962)), fn(963)), fn(971)), fn(982)), fn(195)), fn(972)), fn(973)), fn(974)), fn(975)), fn(978)), fn(979)), fn(980)), fn(981)), fn(983)), fn(984)), fn(995)), fn(985)), fn(986)), fn(988)), fn(989)), fn(991)), fn(992)), fn(993)), fn(994)), fn(996)), fn(997)), fn(998)), fn(999)), fn(1000)), fn(1002)), fn(1003)), fn(1004)), fn(872)), fn(1005)), fn(1001)), fn(1006)), fn(1007)), fn(1737)), min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(fn(1361), fn(1091)), fn(957)), fn(958)), fn(959)), fn(960)), fn(961)), fn(962)), fn(963)), fn(971)), fn(982)), fn(195)), fn(972)), fn(973)), fn(974)), fn(975)), fn(978)), fn(979)), fn(980)), fn(981)), fn(983)), fn(984)), fn(995)), fn(985)), fn(986)), fn(988)), fn(989)), fn(991)), fn(992)), fn(993)), fn(994)), fn(996)), fn(997)), fn(998)), fn(999)), fn(1000)), fn(1002)), fn(1003)), fn(1004)), fn(872)), fn(1005)), fn(1001)), fn(1006)), fn(1008)), fn(1738))), min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(fn(1361), fn(1091)), fn(957)), fn(958)), fn(959)), fn(960)), fn(961)), fn(962)), fn(963)), fn(971)), fn(982)), fn(195)), fn(972)), fn(973)), fn(974)), fn(975)), fn(978)), fn(979)), fn(980)), fn(981)), fn(983)), fn(984)), fn(995)), fn(985)), fn(986)), fn(988)), fn(989)), fn(990)), fn(992)), fn(993)), fn(994)), fn(996)), fn(997)), fn(998)), fn(999)), fn(1000)), fn(1002)), fn(1003)), fn(1004)), fn(872)), fn(1005)), fn(1001)), fn(1006)), fn(1007)), fn(1737))), min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(fn(1361), fn(1091)), fn(957)), fn(958)), fn(959)), fn(960)), fn(961)), fn(962)), fn(963)), fn(971)), fn(982)), fn(195)), fn(972)), fn(973)), fn(974)), fn(975)), fn(978)), fn(979)), fn(980)), fn(981)), fn(983)), fn(984)), fn(995)), fn(985)), fn(986)), fn(988)), fn(989)), fn(990)), fn(992)), fn(993)), fn(994)), fn(996)), fn(997)), fn(998)), fn(999)), fn(1000)), fn(1002)), fn(1003)), fn(1004)), fn(872)), fn(1005)), fn(1001)), fn(1006)), fn(1008)), fn(1737))), min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(fn(1361), fn(1091)), fn(957)), fn(958)), fn(959)), fn(960)), fn(961)), fn(962)), fn(963)), fn(971)), fn(982)), fn(195)), fn(972)), fn(973)), fn(974)), fn(975)), fn(978)), fn(979)), fn(980)), fn(981)), fn(983)), fn(984)), fn(995)), fn(985)), fn(986)), fn(988)), fn(989)), fn(990)), fn(992)), fn(993)), fn(994)), fn(996)), fn(997)), fn(998)), fn(999)), fn(1000)), fn(1002)), fn(1003)), fn(1004)), fn(872)), fn(1005)), fn(1001)), fn(1006)), fn(1007)), fn(1738))), min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(fn(1361), fn(1091)), fn(957)), fn(958)), fn(959)), fn(960)), fn(961)), fn(962)), fn(963)), fn(971)), fn(982)), fn(195)), fn(972)), fn(973)), fn(974)), fn(975)), fn(978)), fn(979)), fn(980)), fn(981)), fn(983)), fn(984)), fn(995)), fn(985)), fn(986)), fn(988)), fn(989)), fn(991)), fn(992)), fn(993)), fn(994)), fn(996)), fn(997)), fn(998)), fn(999)), fn(1000)), fn(1002)), fn(1003)), fn(1004)), fn(872)), fn(1005)), fn(1001)), fn(1006)), fn(1008)), fn(1737))), min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(fn(1361), fn(1091)), fn(957)), fn(958)), fn(959)), fn(960)), fn(961)), fn(962)), fn(963)), fn(971)), fn(982)), fn(195)), fn(972)), fn(973)), fn(974)), fn(975)), fn(978)), fn(979)), fn(980)), fn(981)), fn(983)), fn(984)), fn(995)), fn(985)), fn(986)), fn(988)), fn(989)), fn(990)), fn(992)), fn(993)), fn(994)), fn(996)), fn(997)), fn(998)), fn(999)), fn(1000)), fn(1002)), fn(1003)), fn(1004)), fn(872)), fn(1005)), fn(1001)), fn(1006)), fn(1008)), fn(1738))), min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(fn(1361), fn(1091)), fn(957)), fn(958)), fn(959)), fn(960)), fn(961)), fn(962)), fn(963)), fn(971)), fn(982)), fn(195)), fn(972)), fn(973)), fn(974)), fn(975)), fn(978)), fn(979)), fn(980)), fn(981)), fn(983)), fn(984)), fn(995)), fn(985)), fn(986)), fn(988)), fn(989)), fn(991)), fn(992)), fn(993)), fn(994)), fn(996)), fn(997)), fn(998)), fn(999)), fn(1000)), fn(1002)), fn(1003)), fn(1004)), fn(872)), fn(1005)), fn(1001)), fn(1006)), fn(1007)), fn(1738)))
	v[367] = max(min(min(min(min(min(min(min(min(min(fn(261), fn(954)), fn(739)), fn(648)), fn(94)), fn(1688)), fn(1689)), fn(1690)), fn(1691)), fn(1692)), min(min(min(min(min(min(min(min(min(fn(261), fn(954)), fn(738)), fn(648)), fn(94)), fn(1688)), fn(1689)), fn(1690)), fn(1691)), fn(1692)))
	v[368] = max(min(min(min(min(min(min(min(min(min(fn(1689), fn(1692)), fn(954)), fn(94)), fn(739)), fn(1690)), fn(261)), fn(1691)), fn(648)), fn(1688)), min(min(min(min(min(min(min(min(min(fn(1689), fn(1692)), fn(738)), fn(94)), fn(954)), fn(1690)), fn(1691)), fn(261)), fn(648)), fn(1688)))
	v[369] = default()
	v[370] = max(fn(1481), fn(1918))
	v[371] = default()
	v[372] = default()
	v[373] = default()
	v[374] = default()
	v[375] = default()
	v[376] = default()
	v[377] = default()
	v[378] = default()
	v[379] = default()
	
	# I had to add these ... 
	v[380] = default()
	v[381] = default()
	v[382] = default()
	v[383] = default()
	v[384] = default()
	v[385] = default()
	v[386] = default()
	v[387] = default()
	v[388] = default()

	# growth -
	v[389] = 1.0	# growth -

	return v
end
