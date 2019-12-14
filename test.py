import networkx as nx
import time
import ELOD


def show(G, subgraph=None):
    import matplotlib.pyplot as plt
    nx.draw(G, with_labels=True, node_size=150, width=1.0,
            node_color=[1 if subgraph is None or v in subgraph else 0 for v in G.nodes()],
            edge_color=[1 if subgraph is None or subgraph.has_edge(v,u) else 0 for v,u in G.edges()],
            node_cmap=plt.cm.rainbow,
            edge_cmap=plt.cm.rainbow )
    plt.show()


#evaluate GW https://arxiv.org/pdf/1908.05767.pdf, optimal value 0.7632


core_trace_conductance = list()
core_pcst_conductance = list()
core_trace_a = list()
core_pcst_a = list()
core_trace_time = list()
core_pcst_time = list()
trace_time = list()
pcst_time = list()
trace_ELOD = list()
pcst_ELOD = list()


for i in range(100):
    G = nx.erdos_renyi_graph(100, 0.1, seed=i).to_directed()
    G.remove_nodes_from(list(nx.isolates(G)))
    r = list(G.nodes)[0]

    tic = time.time()
    T, trace_a = ELOD.core(G, r, ELOD.trace)
    core_trace_time.append(time.time()-tic)
    tic = time.time()
    W, pcst_a = ELOD.core(G, r, ELOD.pcst)
    core_pcst_time.append(time.time()-tic)

    core_trace_conductance.append(ELOD.conductance(G, T))
    core_pcst_conductance.append(ELOD.conductance(G, W))

    core_trace_a.append(trace_a)
    core_pcst_a.append(pcst_a)

    found_a = min(trace_a, pcst_a)

    tic = time.time()
    T = ELOD.trace(G, r, found_a)
    trace_time.append(time.time()-tic)
    tic = time.time()
    W = ELOD.pcst(G, r, found_a)
    pcst_time.append(time.time()-tic)

    trace_ELOD.append(ELOD.ELOD(G, T, found_a))
    pcst_ELOD.append(ELOD.ELOD(G, W, found_a))

    print(i)

print('Core trace conductance')
print(core_trace_conductance)
print('Core PCST conductance')
print(core_pcst_conductance)

print('Core trace a')
print(core_trace_a)
print('Core PCST a')
print(core_pcst_a)

print('Core trace time')
print(core_trace_time)
print('Core PCST time')
print(core_trace_time)

print('Trace time')
print(trace_time)
print('PCST time')
print(pcst_time)

print('Trace ELOD')
print(trace_ELOD)
print('PCST ELOD')
print(pcst_ELOD)

"""
EXAMPLE RESULS (SIZE=100)
Core trace conductance
[17.0, 15.7, 15.833333333333332, 18.5, 19.5, 14.88888888888889, 19.666666666666668, 14.5, 19.333333333333332, 17.75, 17.5, 30.0, 14.5, 15.75, 16.285714285714285, 16.166666666666668, 28.0, 13.5, 17.5, 17.2, 19.25, 14.363636363636363, 18.75, 16.285714285714285, 13.0, 21.5, 18.5, 17.0, 15.777777777777779, 18.0, 14.0, 13.5, 16.333333333333332, 16.666666666666668, 15.5, 16.75, 18.4, 18.0, 16.666666666666668, 15.333333333333332, 17.75, 17.0, 17.0, 30.0, 19.5, 33.0, 17.0, 16.0, 15.399999999999999, 19.5, 14.88888888888889, 15.5, 14.666666666666666, 15.5, 18.666666666666668, 16.75, 17.8, 21.333333333333332, 15.166666666666668, 13.833333333333334, 24.0, 18.0, 15.0, 17.0, 17.6, 15.571428571428573, 14.384615384615385, 18.8, 20.0, 29.0, 14.454545454545455, 16.4, 20.5, 15.222222222222221, 15.0, 17.0, 15.833333333333332, 18.4, 21.0, 15.571428571428573, 18.25, 15.166666666666668, 16.625, 17.5, 15.875, 15.666666666666668, 30.0, 17.25, 15.666666666666668, 16.166666666666668, 15.625, 16.0, 16.333333333333332, 19.4, 14.875, 14.777777777777779, 16.0, 16.142857142857142, 16.5, 15.0]
Core PCST conductance
[16.0, 13.714285714285714, 11.692307692307692, 12.923076923076923, 18.5, 14.285714285714286, 12.11111111111111, 14.0, 14.555555555555555, 14.066666666666666, 12.352941176470589, 12.538461538461538, 17.666666666666668, 12.538461538461538, 11.75609756097561, 12.428571428571429, 11.681818181818182, 10.535714285714286, 18.0, 13.23076923076923, 14.0, 12.473684210526315, 13.3, 14.0, 10.181818181818182, 13.176470588235293, 15.100000000000001, 15.399999999999999, 15.285714285714285, 13.7, 11.466666666666667, 12.9, 13.9, 14.142857142857142, 12.466666666666667, 13.933333333333334, 17.2, 13.444444444444445, 14.076923076923077, 11.9, 12.5, 13.529411764705882, 12.777777777777779, 12.0, 12.869565217391305, 11.64516129032258, 13.5, 12.5, 14.25, 16.666666666666668, 13.222222222222221, 12.2, 12.115384615384615, 12.294117647058824, 17.0, 13.0, 27.0, 14.11111111111111, 12.23076923076923, 12.333333333333334, 23.0, 13.714285714285714, 21.0, 13.0, 17.0, 12.777777777777779, 14.8, 28.0, 15.666666666666668, 11.72, 13.8, 10.3125, 15.333333333333332, 12.142857142857142, 12.416666666666666, 12.90625, 11.689655172413794, 12.882352941176471, 15.833333333333332, 24.0, 13.272727272727273, 14.666666666666666, 14.153846153846153, 12.38888888888889, 12.1, 14.6, 11.772727272727273, 10.838709677419354, 11.807692307692308, 12.0, 12.23076923076923, 13.833333333333334, 14.0, 18.2, 13.875, 20.0, 15.25, 26.0, 13.357142857142858, 11.555555555555555]
Core trace a
[16.714285714285715, 14.576923076923077, 14.31578947368421, 14.238095238095237, 15.142857142857142, 13.782608695652174, 14.666666666666666, 14.0, 16.555555555555557, 16.470588235294116, 13.736842105263158, 17.0, 13.217391304347826, 14.833333333333334, 14.380952380952381, 14.333333333333334, 16.0, 12.761904761904763, 16.0, 14.882352941176471, 15.866666666666667, 13.833333333333334, 15.117647058823529, 14.88888888888889, 12.74074074074074, 17.0, 15.333333333333334, 14.818181818181818, 14.347826086956522, 16.625, 12.956521739130435, 13.0, 13.857142857142858, 14.882352941176471, 14.857142857142858, 14.84, 15.588235294117647, 13.789473684210526, 14.80952380952381, 14.23076923076923, 14.705882352941176, 15.157894736842104, 14.333333333333334, 17.571428571428573, 15.3, 15.0, 14.555555555555555, 14.588235294117647, 13.11111111111111, 16.22222222222222, 14.380952380952381, 14.380952380952381, 13.695652173913043, 14.238095238095237, 15.285714285714286, 15.5, 15.428571428571429, 16.533333333333335, 14.0, 12.272727272727273, 15.181818181818182, 16.857142857142858, 13.571428571428571, 14.210526315789474, 15.105263157894736, 14.333333333333334, 13.96, 15.7, 16.181818181818183, 16.22222222222222, 13.379310344827585, 14.214285714285714, 16.142857142857142, 14.391304347826088, 13.55, 15.130434782608695, 13.857142857142858, 14.65, 16.0, 14.0, 15.928571428571429, 13.52, 15.933333333333334, 14.375, 15.0, 15.0625, 17.166666666666668, 13.04, 14.277777777777779, 13.6, 14.842105263157896, 14.6875, 15.0, 15.761904761904763, 14.117647058823529, 13.96, 13.047619047619047, 14.647058823529411, 14.875, 13.08695652173913]
Core PCST a
[13.75, 12.63265306122449, 11.76, 12.448275862068966, 13.25, 12.166666666666666, 12.264150943396226, 12.423076923076923, 12.535714285714286, 13.018518518518519, 11.422222222222222, 12.228070175438596, 12.413793103448276, 12.326530612244898, 11.830508474576272, 11.285714285714286, 11.470588235294118, 11.16326530612245, 14.0, 12.724137931034482, 14.166666666666666, 12.0, 12.382978723404255, 12.72, 9.353535353535353, 13.428571428571429, 11.89090909090909, 13.727272727272727, 14.266666666666667, 13.346153846153847, 11.352941176470589, 12.73913043478261, 12.73076923076923, 13.12, 12.081632653061224, 12.808510638297872, 15.846153846153847, 12.478260869565217, 12.787234042553191, 12.25531914893617, 12.0, 12.605263157894736, 11.471698113207546, 12.196078431372548, 12.11864406779661, 11.457627118644067, 12.586206896551724, 12.758620689655173, 12.142857142857142, 13.545454545454545, 12.588235294117647, 11.942307692307692, 11.826923076923077, 11.6, 14.153846153846153, 12.280701754385966, 15.142857142857142, 12.714285714285714, 12.37037037037037, 11.241379310344827, 14.222222222222221, 13.625, 11.217391304347826, 12.0, 13.275862068965518, 12.580645161290322, 13.037037037037036, 14.714285714285714, 12.857142857142858, 12.085714285714285, 12.393939393939394, 10.836363636363636, 12.380952380952381, 12.294117647058824, 12.214285714285714, 12.627118644067796, 11.833333333333334, 12.057692307692308, 11.6875, 15.333333333333334, 12.478260869565217, 13.15, 12.653061224489797, 12.104166666666666, 11.568627450980392, 13.217391304347826, 11.14516129032258, 11.218181818181819, 11.771929824561404, 11.475409836065573, 12.488372093023257, 13.19047619047619, 12.428571428571429, 14.68, 13.107142857142858, 13.923076923076923, 12.9, 13.9, 12.6, 11.583333333333334]
Core trace time
[0.09719109535217285, 0.08182096481323242, 0.10735845565795898, 0.08004093170166016, 0.10745382308959961, 0.07758259773254395, 0.09709930419921875, 0.08398294448852539, 0.1156461238861084, 0.08597588539123535, 0.11735796928405762, 0.10374164581298828, 0.10153722763061523, 0.09726500511169434, 0.1055452823638916, 0.10272669792175293, 0.13098359107971191, 0.07735562324523926, 0.11211228370666504, 0.0957944393157959, 0.10233616828918457, 0.09309673309326172, 0.1377544403076172, 0.08410191535949707, 0.11506271362304688, 0.10069775581359863, 0.08108830451965332, 0.11412620544433594, 0.08284139633178711, 0.10260295867919922, 0.0819551944732666, 0.1117093563079834, 0.0828549861907959, 0.09768295288085938, 0.08179783821105957, 0.10433697700500488, 0.09546661376953125, 0.12063193321228027, 0.09618663787841797, 0.10701584815979004, 0.12773394584655762, 0.13750648498535156, 0.08478474617004395, 0.11896228790283203, 0.11913919448852539, 0.08315396308898926, 0.10249209403991699, 0.08086490631103516, 0.1093297004699707, 0.11472129821777344, 0.11268258094787598, 0.08219385147094727, 0.10950469970703125, 0.08517980575561523, 0.11458063125610352, 0.11880660057067871, 0.10139632225036621, 0.08292484283447266, 0.09394049644470215, 0.10434675216674805, 0.09494519233703613, 0.12201809883117676, 0.07695794105529785, 0.11122369766235352, 0.09791016578674316, 0.09839653968811035, 0.08577799797058105, 0.09654641151428223, 0.09914493560791016, 0.11000800132751465, 0.0946967601776123, 0.09044122695922852, 0.0915372371673584, 0.11290526390075684, 0.0776820182800293, 0.10072898864746094, 0.08128142356872559, 0.0958409309387207, 0.09595203399658203, 0.09695553779602051, 0.09625101089477539, 0.11124682426452637, 0.08280396461486816, 0.09596610069274902, 0.08072233200073242, 0.0937349796295166, 0.09755373001098633, 0.09215259552001953, 0.08158111572265625, 0.10994982719421387, 0.07992887496948242, 0.09635615348815918, 0.08068704605102539, 0.10010361671447754, 0.08054995536804199, 0.09432792663574219, 0.07656383514404297, 0.1081228256225586, 0.07895445823669434, 0.09360599517822266]
Core PCST time
[0.09719109535217285, 0.08182096481323242, 0.10735845565795898, 0.08004093170166016, 0.10745382308959961, 0.07758259773254395, 0.09709930419921875, 0.08398294448852539, 0.1156461238861084, 0.08597588539123535, 0.11735796928405762, 0.10374164581298828, 0.10153722763061523, 0.09726500511169434, 0.1055452823638916, 0.10272669792175293, 0.13098359107971191, 0.07735562324523926, 0.11211228370666504, 0.0957944393157959, 0.10233616828918457, 0.09309673309326172, 0.1377544403076172, 0.08410191535949707, 0.11506271362304688, 0.10069775581359863, 0.08108830451965332, 0.11412620544433594, 0.08284139633178711, 0.10260295867919922, 0.0819551944732666, 0.1117093563079834, 0.0828549861907959, 0.09768295288085938, 0.08179783821105957, 0.10433697700500488, 0.09546661376953125, 0.12063193321228027, 0.09618663787841797, 0.10701584815979004, 0.12773394584655762, 0.13750648498535156, 0.08478474617004395, 0.11896228790283203, 0.11913919448852539, 0.08315396308898926, 0.10249209403991699, 0.08086490631103516, 0.1093297004699707, 0.11472129821777344, 0.11268258094787598, 0.08219385147094727, 0.10950469970703125, 0.08517980575561523, 0.11458063125610352, 0.11880660057067871, 0.10139632225036621, 0.08292484283447266, 0.09394049644470215, 0.10434675216674805, 0.09494519233703613, 0.12201809883117676, 0.07695794105529785, 0.11122369766235352, 0.09791016578674316, 0.09839653968811035, 0.08577799797058105, 0.09654641151428223, 0.09914493560791016, 0.11000800132751465, 0.0946967601776123, 0.09044122695922852, 0.0915372371673584, 0.11290526390075684, 0.0776820182800293, 0.10072898864746094, 0.08128142356872559, 0.0958409309387207, 0.09595203399658203, 0.09695553779602051, 0.09625101089477539, 0.11124682426452637, 0.08280396461486816, 0.09596610069274902, 0.08072233200073242, 0.0937349796295166, 0.09755373001098633, 0.09215259552001953, 0.08158111572265625, 0.10994982719421387, 0.07992887496948242, 0.09635615348815918, 0.08068704605102539, 0.10010361671447754, 0.08054995536804199, 0.09432792663574219, 0.07656383514404297, 0.1081228256225586, 0.07895445823669434, 0.09360599517822266]
Trace time
[0.015499591827392578, 0.015778541564941406, 0.015150308609008789, 0.015233278274536133, 0.015450716018676758, 0.015015840530395508, 0.015431880950927734, 0.017686843872070312, 0.015199661254882812, 0.016122102737426758, 0.017668962478637695, 0.0168917179107666, 0.014128446578979492, 0.015479326248168945, 0.016090869903564453, 0.015832901000976562, 0.017169952392578125, 0.014871835708618164, 0.015388965606689453, 0.015473127365112305, 0.015690326690673828, 0.016423702239990234, 0.015538692474365234, 0.015786170959472656, 0.015029191970825195, 0.016205310821533203, 0.016339778900146484, 0.029709339141845703, 0.020290613174438477, 0.01618194580078125, 0.01598072052001953, 0.015306711196899414, 0.015078544616699219, 0.02127242088317871, 0.015732288360595703, 0.01795792579650879, 0.02862715721130371, 0.01423501968383789, 0.015869855880737305, 0.027926921844482422, 0.016597509384155273, 0.015822410583496094, 0.0167696475982666, 0.02207326889038086, 0.026139020919799805, 0.016479969024658203, 0.016508817672729492, 0.015496253967285156, 0.014465570449829102, 0.01742696762084961, 0.016328096389770508, 0.01850104331970215, 0.015169620513916016, 0.01761460304260254, 0.016272544860839844, 0.016227245330810547, 0.015905380249023438, 0.01717829704284668, 0.014971256256103516, 0.014038324356079102, 0.01510167121887207, 0.014976978302001953, 0.014712810516357422, 0.01521754264831543, 0.01557016372680664, 0.01602935791015625, 0.015115976333618164, 0.029300212860107422, 0.01616954803466797, 0.015211343765258789, 0.015049934387207031, 0.015333414077758789, 0.01455831527709961, 0.016163349151611328, 0.014535188674926758, 0.016384601593017578, 0.015472412109375, 0.014947891235351562, 0.01508188247680664, 0.014777660369873047, 0.01593470573425293, 0.015040159225463867, 0.016003131866455078, 0.02881646156311035, 0.014615058898925781, 0.015137434005737305, 0.01502084732055664, 0.01559305191040039, 0.01535487174987793, 0.01581096649169922, 0.014989614486694336, 0.01579451560974121, 0.014687776565551758, 0.016678810119628906, 0.01562952995300293, 0.015590429306030273, 0.014671802520751953, 0.015290498733520508, 0.014609813690185547, 0.01483464241027832]
PCST time
[0.008530139923095703, 0.008467674255371094, 0.008152484893798828, 0.008152246475219727, 0.008132457733154297, 0.008200407028198242, 0.008763313293457031, 0.014584541320800781, 0.008002996444702148, 0.013650178909301758, 0.009737014770507812, 0.008646726608276367, 0.007848024368286133, 0.008337259292602539, 0.009561538696289062, 0.009305477142333984, 0.008622884750366211, 0.008213281631469727, 0.008315801620483398, 0.00836801528930664, 0.008814334869384766, 0.009354114532470703, 0.008667230606079102, 0.008962869644165039, 0.008434057235717773, 0.008743762969970703, 0.008440017700195312, 0.009856700897216797, 0.00965428352355957, 0.009505748748779297, 0.010543107986450195, 0.008153676986694336, 0.009005069732666016, 0.008640289306640625, 0.008892536163330078, 0.0192110538482666, 0.013205528259277344, 0.008638858795166016, 0.008689641952514648, 0.010276556015014648, 0.03345775604248047, 0.008861541748046875, 0.010850906372070312, 0.009685516357421875, 0.017264604568481445, 0.008859395980834961, 0.009332656860351562, 0.008261680603027344, 0.008201360702514648, 0.00911712646484375, 0.014230966567993164, 0.010849475860595703, 0.008255481719970703, 0.008514165878295898, 0.008711576461791992, 0.009032726287841797, 0.00931859016418457, 0.017298221588134766, 0.007842302322387695, 0.007453203201293945, 0.008373737335205078, 0.008112430572509766, 0.007582426071166992, 0.008480072021484375, 0.008882284164428711, 0.008868694305419922, 0.00830221176147461, 0.009228229522705078, 0.009217023849487305, 0.00831747055053711, 0.008184194564819336, 0.008523941040039062, 0.007904291152954102, 0.008665323257446289, 0.007989883422851562, 0.010271310806274414, 0.008390426635742188, 0.008331298828125, 0.008205652236938477, 0.0075283050537109375, 0.008283376693725586, 0.008419036865234375, 0.008597373962402344, 0.008480310440063477, 0.007842540740966797, 0.008308649063110352, 0.008791923522949219, 0.008758544921875, 0.008817434310913086, 0.011794567108154297, 0.00876927375793457, 0.00821828842163086, 0.008275270462036133, 0.008572101593017578, 0.008130311965942383, 0.008322954177856445, 0.007716655731201172, 0.007843494415283203, 0.008182287216186523, 0.008289337158203125]
Trace ELOD
[20.25, 50.551020408163254, 48.0, 37.58620689655169, 28.5, 39.0, 43.24528301886792, 34.69230769230768, 37.571428571428555, 60.611111111111086, 43.97777777777776, 57.298245614035125, 25.448275862068954, 32.83673469387753, 56.59322033898303, 48.28571428571428, 38.05882352941177, 33.571428571428555, 12.0, 36.68965517241381, 32.83333333333334, 44.0, 46.72340425531914, 39.03999999999999, 121.26262626262627, 28.571428571428555, 70.72727272727275, 31.454545454545467, 22.599999999999994, 29.5, 36.882352941176464, 21.13043478260869, 30.769230769230774, 31.680000000000007, 38.36734693877554, 50.78723404255322, 19.058823529411768, 28.739130434782624, 42.468085106383, 51.361702127659555, 46.0, 49.10526315789474, 46.84905660377359, 53.49019607843138, 62.27118644067798, 55.72881355932202, 35.79310344827587, 31.103448275862064, 19.571428571428584, 25.0, 37.64705882352939, 51.730769230769226, 42.980769230769226, 58.19999999999999, 20.769230769230774, 57.859649122807014, 18.285714285714292, 57.0, 28.18518518518519, 26.103448275862064, 12.111111111111114, 23.5, 48.65217391304347, 42.0, 36.13793103448276, 36.80645161290323, 30.518518518518533, 27.714285714285722, 37.714285714285694, 45.94285714285718, 36.333333333333314, 51.74545454545455, 28.19047619047619, 46.64705882352939, 27.571428571428584, 52.4406779661017, 42.5, 51.96153846153845, 50.5625, 18.0, 50.0, 19.5, 51.63265306122446, 36.33333333333334, 48.0392156862745, 30.17391304347825, 60.935483870967744, 45.545454545454504, 46.017543859649095, 51.63934426229508, 44.72093023255812, 24.904761904761926, 35.28571428571428, 27.19999999999999, 24.821428571428555, 16.692307692307693, 16.89999999999999, 23.0, 36.599999999999994, 34.91666666666666]
PCST ELOD
[6.5, 29.14285714285714, 24.24000000000001, 19.17241379310343, 12.5, 21.833333333333343, 22.8679245283019, 20.615384615384613, 27.178571428571416, 30.72222222222223, 32.82222222222222, 34.07017543859649, 18.758620689655174, 31.51020408163265, 37.94915254237287, 15.0, 26.64705882352942, 10.428571428571445, 5.0, 19.586206896551744, 14.166666666666686, 28.0, 38.34042553191489, 13.679999999999993, 100.55555555555554, 12.714285714285722, 42.09090909090909, 13.363636363636374, 14.133333333333326, 13.538461538461547, 16.70588235294116, 11.608695652173907, 21.692307692307708, 14.16000000000001, 20.775510204081655, 31.872340425531917, 14.058823529411768, 17.695652173913047, 29.765957446808528, 19.34042553191489, 15.0, 32.71052631578948, 20.754716981132077, 20.901960784313758, 40.27118644067798, 36.81355932203388, 22.965517241379303, 10.379310344827587, 12.42857142857143, 12.363636363636367, 14.705882352941174, 31.442307692307736, 33.5, 28.80000000000001, 15.384615384615387, 34.385964912280684, 12.857142857142858, 43.14285714285717, 11.18518518518519, 12.551724137931032, 9.777777777777779, 7.625, 10.782608695652174, 14.0, 14.172413793103445, 10.774193548387103, 13.814814814814824, 14.285714285714286, 11.42857142857143, 15.85714285714289, 12.030303030303031, 15.236363636363649, 11.857142857142854, 23.764705882352928, 14.428571428571445, 40.93220338983053, 24.833333333333314, 31.019230769230745, 30.875, 0.0, 39.47826086956525, 15.099999999999994, 32.51020408163265, 23.125, 30.627450980392155, 11.913043478260875, 35.80645161290323, 19.23636363636365, 26.929824561403507, 36.59016393442624, 19.302325581395337, 9.857142857142861, 10.285714285714285, 22.599999999999994, 14.142857142857139, 7.076923076923077, 13.399999999999999, 13.1, 24.599999999999994, 17.5]


EXAMPLE RESULS (SIZE=200)
Core trace conductance
[29.857142857142858, 38.333333333333336, 39.0, 32.25, 26.0, 30.428571428571427, 34.8, 35.333333333333336, 41.0, 48.0, 39.0, 55.0, 52.0, 38.0, 29.857142857142858, 29.375, 30.333333333333332, 33.8, 31.125, 33.8, 28.4, 38.0, 57.0, 33.6, 36.333333333333336, 34.0, 28.666666666666668, 32.0, 32.2, 31.333333333333336, 29.571428571428573, 31.166666666666664, 30.818181818181817, 35.0, 31.166666666666664, 32.714285714285715, 29.2, 60.0, 29.666666666666668, 26.166666666666668, 31.444444444444443, 31.5, 33.75, 31.625, 34.5, 58.0, 39.666666666666664, 33.2, 40.0, 28.894736842105264, 53.0, 31.0, 30.5, 35.5, 29.833333333333332, 34.25, 31.6, 32.8, 29.857142857142858, 33.5, 31.0, 29.5, 29.857142857142858, 34.25, 38.0, 30.666666666666668, 26.529411764705884, 36.0, 31.799999999999997, 29.22222222222222, 58.0, 30.625, 31.333333333333336, 32.4, 31.5, 33.5, 41.5, 48.5, 30.8, 28.75, 35.666666666666664, 31.799999999999997, 34.2, 27.285714285714285, 43.0, 39.5, 26.764705882352942, 33.25, 30.428571428571427, 28.4375, 34.666666666666664, 28.857142857142858, 42.0, 34.0, 31.0, 36.666666666666664, 27.77777777777778, 29.0, 40.0, 28.428571428571427]
Core PCST conductance
[26.833333333333332, 24.666666666666668, 34.0, 25.11111111111111, 24.571428571428573, 23.361702127659573, 50.0, 22.46153846153846, 24.58974358974359, 59.0, 27.22222222222222, 31.799999999999997, 24.15625, 47.0, 24.77777777777778, 24.68, 23.236363636363638, 25.25, 23.428571428571427, 25.5625, 23.5, 28.0, 24.25, 27.53846153846154, 23.389830508474578, 24.523809523809526, 24.488888888888887, 24.26829268292683, 24.405405405405407, 26.5, 25.4, 24.193548387096776, 27.454545454545453, 28.833333333333332, 25.571428571428573, 25.94736842105263, 26.22222222222222, 22.566037735849058, 47.0, 23.4, 26.6, 24.857142857142858, 50.0, 24.142857142857142, 23.942857142857143, 24.32, 23.476190476190474, 26.357142857142858, 24.0, 28.5, 22.076923076923077, 24.636363636363637, 25.470588235294116, 23.733333333333334, 26.666666666666668, 27.333333333333332, 25.58823529411765, 23.583333333333332, 24.736842105263158, 24.96551724137931, 28.125, 44.0, 24.62962962962963, 25.59090909090909, 51.0, 26.0, 23.944444444444443, 28.416666666666668, 24.904761904761905, 35.0, 24.685714285714287, 25.6, 25.72222222222222, 23.586206896551722, 25.307692307692307, 25.0, 34.333333333333336, 25.916666666666668, 26.571428571428573, 22.416666666666668, 24.607142857142858, 28.285714285714285, 26.208333333333332, 24.0, 26.0, 26.4, 27.166666666666668, 27.428571428571427, 27.833333333333332, 25.25925925925926, 23.925, 28.0, 25.0, 24.68, 27.75, 24.94736842105263, 27.714285714285715, 25.636363636363637, 22.145454545454545, 28.4]
Core trace a
[28.6, 29.88235294117647, 28.352941176470587, 27.41176470588235, 25.0, 28.333333333333332, 30.083333333333332, 27.35, 31.076923076923077, 32.55555555555556, 29.11111111111111, 28.136363636363637, 28.923076923076923, 29.266666666666666, 27.705882352941178, 27.037037037037038, 28.625, 29.142857142857142, 28.88235294117647, 29.866666666666667, 26.692307692307693, 29.0625, 36.25, 30.285714285714285, 30.5, 28.894736842105264, 26.38095238095238, 27.944444444444443, 28.764705882352942, 28.38888888888889, 27.761904761904763, 28.45, 27.235294117647058, 29.666666666666668, 28.055555555555557, 29.176470588235293, 27.863636363636363, 29.857142857142858, 26.7027027027027, 24.695652173913043, 27.526315789473685, 26.666666666666668, 27.954545454545453, 29.72222222222222, 28.928571428571427, 31.3, 30.11764705882353, 28.454545454545453, 27.35, 26.94736842105263, 28.133333333333333, 27.227272727272727, 27.666666666666668, 27.8125, 26.708333333333332, 30.0, 28.19047619047619, 29.928571428571427, 27.77777777777778, 28.181818181818183, 29.0, 27.227272727272727, 28.0, 29.666666666666668, 30.5, 28.58823529411765, 25.130434782608695, 29.58823529411765, 28.529411764705884, 27.45, 39.666666666666664, 27.95, 28.333333333333332, 27.68421052631579, 27.647058823529413, 30.055555555555557, 29.307692307692307, 31.77777777777778, 27.57894736842105, 25.023255813953487, 29.307692307692307, 28.928571428571427, 28.608695652173914, 25.8, 30.46153846153846, 28.045454545454547, 25.5, 28.5, 28.523809523809526, 26.024390243902438, 29.428571428571427, 27.318181818181817, 30.066666666666666, 28.45, 27.11764705882353, 29.642857142857142, 26.94736842105263, 26.73913043478261, 27.263157894736842, 26.333333333333332]
Core PCST a
[25.208333333333332, 23.537735849056602, 23.112903225806452, 23.83050847457627, 22.76923076923077, 22.44144144144144, 25.61111111111111, 22.39047619047619, 23.285714285714285, 26.217391304347824, 24.979166666666668, 26.041666666666668, 24.48076923076923, 27.8, 23.74468085106383, 24.5, 23.19191919191919, 24.326923076923077, 23.019417475728154, 24.73913043478261, 23.083333333333332, 23.8, 24.2, 25.346938775510203, 23.31, 24.220338983050848, 23.33, 22.63, 24.535714285714285, 25.181818181818183, 24.770833333333332, 24.416666666666668, 24.98181818181818, 24.979166666666668, 24.036363636363635, 25.568181818181817, 25.333333333333332, 22.285714285714285, 25.285714285714285, 21.860869565217392, 25.306122448979593, 24.169811320754718, 24.71153846153846, 23.73404255319149, 22.546296296296298, 23.288659793814432, 22.88888888888889, 25.0, 22.484848484848484, 25.07843137254902, 21.96330275229358, 24.181818181818183, 24.51851851851852, 23.16923076923077, 23.75, 25.06, 24.925, 22.635514018691588, 24.434782608695652, 24.9375, 26.318181818181817, 29.333333333333332, 24.076923076923077, 25.666666666666668, 27.384615384615383, 24.24074074074074, 22.694444444444443, 25.318181818181817, 24.642857142857142, 23.9375, 23.86046511627907, 24.555555555555557, 25.03921568627451, 23.842105263157894, 23.686274509803923, 23.99056603773585, 27.0, 24.03846153846154, 24.304347826086957, 22.11111111111111, 23.26530612244898, 24.06153846153846, 25.666666666666668, 23.714285714285715, 23.925925925925927, 25.133333333333333, 24.44, 23.73076923076923, 25.137931034482758, 25.21951219512195, 23.20618556701031, 25.61904761904762, 24.56, 24.717391304347824, 26.289473684210527, 24.61111111111111, 25.85, 24.90909090909091, 22.042735042735043, 25.157894736842106]
Core trace time
[0.33826375007629395, 0.39800333976745605, 0.35762619972229004, 0.4207417964935303, 0.2899346351623535, 0.38587093353271484, 0.340714693069458, 0.35191988945007324, 0.37697339057922363, 0.37348484992980957, 0.38202667236328125, 0.3657341003417969, 0.3735530376434326, 0.3439016342163086, 0.36620259284973145, 0.3686227798461914, 0.36935853958129883, 0.35222911834716797, 0.3568086624145508, 0.34304118156433105, 0.3059709072113037, 0.3515293598175049, 0.41495752334594727, 0.35285043716430664, 0.359205961227417, 0.3599574565887451, 0.2996561527252197, 0.34845447540283203, 0.36115360260009766, 0.3791680335998535, 0.3584299087524414, 0.3688333034515381, 0.3047609329223633, 0.3642423152923584, 0.36732053756713867, 0.3482937812805176, 0.3732016086578369, 0.347881555557251, 0.3090987205505371, 0.3029477596282959, 0.3250393867492676, 0.3628714084625244, 0.3692154884338379, 0.35514235496520996, 0.3678860664367676, 0.35663795471191406, 0.36822056770324707, 0.37961435317993164, 0.35692620277404785, 0.3025496006011963, 0.3574857711791992, 0.36686110496520996, 0.3670384883880615, 0.34916090965270996, 0.3577699661254883, 0.3527092933654785, 0.35642194747924805, 0.3639066219329834, 0.3710951805114746, 0.3545830249786377, 0.3623192310333252, 0.36064767837524414, 0.3516385555267334, 0.36032867431640625, 0.3569021224975586, 0.36193251609802246, 0.305072546005249, 0.35738563537597656, 0.3553335666656494, 0.35495877265930176, 0.4159703254699707, 0.37552785873413086, 0.3599729537963867, 0.35216403007507324, 0.35692596435546875, 0.3655524253845215, 0.368438720703125, 0.361339807510376, 0.3576054573059082, 0.30227231979370117, 0.364656925201416, 0.3873007297515869, 0.3754153251647949, 0.31459999084472656, 0.3518857955932617, 0.3786029815673828, 0.3216588497161865, 0.3670463562011719, 0.36949777603149414, 0.31599926948547363, 0.3599247932434082, 0.3599874973297119, 0.3577141761779785, 0.3527202606201172, 0.3387575149536133, 0.35162830352783203, 0.34896254539489746, 0.3665015697479248, 0.36922240257263184, 0.3430173397064209]
Core PCST time
[0.33826375007629395, 0.39800333976745605, 0.35762619972229004, 0.4207417964935303, 0.2899346351623535, 0.38587093353271484, 0.340714693069458, 0.35191988945007324, 0.37697339057922363, 0.37348484992980957, 0.38202667236328125, 0.3657341003417969, 0.3735530376434326, 0.3439016342163086, 0.36620259284973145, 0.3686227798461914, 0.36935853958129883, 0.35222911834716797, 0.3568086624145508, 0.34304118156433105, 0.3059709072113037, 0.3515293598175049, 0.41495752334594727, 0.35285043716430664, 0.359205961227417, 0.3599574565887451, 0.2996561527252197, 0.34845447540283203, 0.36115360260009766, 0.3791680335998535, 0.3584299087524414, 0.3688333034515381, 0.3047609329223633, 0.3642423152923584, 0.36732053756713867, 0.3482937812805176, 0.3732016086578369, 0.347881555557251, 0.3090987205505371, 0.3029477596282959, 0.3250393867492676, 0.3628714084625244, 0.3692154884338379, 0.35514235496520996, 0.3678860664367676, 0.35663795471191406, 0.36822056770324707, 0.37961435317993164, 0.35692620277404785, 0.3025496006011963, 0.3574857711791992, 0.36686110496520996, 0.3670384883880615, 0.34916090965270996, 0.3577699661254883, 0.3527092933654785, 0.35642194747924805, 0.3639066219329834, 0.3710951805114746, 0.3545830249786377, 0.3623192310333252, 0.36064767837524414, 0.3516385555267334, 0.36032867431640625, 0.3569021224975586, 0.36193251609802246, 0.305072546005249, 0.35738563537597656, 0.3553335666656494, 0.35495877265930176, 0.4159703254699707, 0.37552785873413086, 0.3599729537963867, 0.35216403007507324, 0.35692596435546875, 0.3655524253845215, 0.368438720703125, 0.361339807510376, 0.3576054573059082, 0.30227231979370117, 0.364656925201416, 0.3873007297515869, 0.3754153251647949, 0.31459999084472656, 0.3518857955932617, 0.3786029815673828, 0.3216588497161865, 0.3670463562011719, 0.36949777603149414, 0.31599926948547363, 0.3599247932434082, 0.3599874973297119, 0.3577141761779785, 0.3527202606201172, 0.3387575149536133, 0.35162830352783203, 0.34896254539489746, 0.3665015697479248, 0.36922240257263184, 0.3430173397064209]
Trace time
[0.0516512393951416, 0.055191755294799805, 0.05385899543762207, 0.06835055351257324, 0.05753684043884277, 0.05423545837402344, 0.056738853454589844, 0.05207681655883789, 0.07226848602294922, 0.053607940673828125, 0.07140970230102539, 0.05258440971374512, 0.06702446937561035, 0.0531463623046875, 0.0515894889831543, 0.0529332160949707, 0.06921577453613281, 0.05178093910217285, 0.05324387550354004, 0.05316305160522461, 0.05266284942626953, 0.05386638641357422, 0.05134177207946777, 0.05427145957946777, 0.06893730163574219, 0.052323102951049805, 0.05175018310546875, 0.05138731002807617, 0.05469870567321777, 0.05400705337524414, 0.05571413040161133, 0.051804542541503906, 0.06350445747375488, 0.05298972129821777, 0.051691532135009766, 0.051310062408447266, 0.05377507209777832, 0.05112600326538086, 0.054575204849243164, 0.052184343338012695, 0.06654596328735352, 0.05312514305114746, 0.05296468734741211, 0.05228447914123535, 0.05180025100708008, 0.054418325424194336, 0.052071571350097656, 0.053119659423828125, 0.07367753982543945, 0.05251121520996094, 0.05144691467285156, 0.05304455757141113, 0.053228139877319336, 0.05198860168457031, 0.05387473106384277, 0.052233219146728516, 0.0527501106262207, 0.05537724494934082, 0.04998159408569336, 0.05190896987915039, 0.052841901779174805, 0.05057644844055176, 0.052987098693847656, 0.06293487548828125, 0.051918983459472656, 0.05422067642211914, 0.05264592170715332, 0.052490949630737305, 0.056725502014160156, 0.05242323875427246, 0.05357670783996582, 0.05072188377380371, 0.05591559410095215, 0.05752825736999512, 0.052445411682128906, 0.05402660369873047, 0.051496267318725586, 0.05072188377380371, 0.05532670021057129, 0.05305218696594238, 0.05315589904785156, 0.07221364974975586, 0.05443859100341797, 0.05197429656982422, 0.05266261100769043, 0.05242443084716797, 0.05119633674621582, 0.05453133583068848, 0.05499672889709473, 0.05281233787536621, 0.052388906478881836, 0.057358503341674805, 0.05399966239929199, 0.052161216735839844, 0.05249595642089844, 0.05176544189453125, 0.05108642578125, 0.04984259605407715, 0.0665125846862793, 0.05342245101928711]
PCST time
[0.03170895576477051, 0.060692548751831055, 0.031177043914794922, 0.048316001892089844, 0.030655622482299805, 0.03390216827392578, 0.028825998306274414, 0.029450654983520508, 0.04826974868774414, 0.029878854751586914, 0.03274273872375488, 0.02884221076965332, 0.043677330017089844, 0.04692888259887695, 0.02755451202392578, 0.029112577438354492, 0.03543496131896973, 0.029407024383544922, 0.029137849807739258, 0.028766393661499023, 0.030671358108520508, 0.03087592124938965, 0.02796339988708496, 0.034857988357543945, 0.030096054077148438, 0.028435945510864258, 0.02878594398498535, 0.04483461380004883, 0.029761791229248047, 0.029178619384765625, 0.030353307723999023, 0.028860807418823242, 0.033927202224731445, 0.029222488403320312, 0.028789997100830078, 0.029459238052368164, 0.030568599700927734, 0.028425931930541992, 0.03035449981689453, 0.028127670288085938, 0.02944350242614746, 0.02884960174560547, 0.028494596481323242, 0.0287020206451416, 0.028840065002441406, 0.03245663642883301, 0.02860260009765625, 0.029541015625, 0.02893352508544922, 0.029635190963745117, 0.02832174301147461, 0.02884364128112793, 0.029047727584838867, 0.03415799140930176, 0.028530359268188477, 0.028748512268066406, 0.030190706253051758, 0.028269052505493164, 0.02849602699279785, 0.029622793197631836, 0.0501251220703125, 0.02796339988708496, 0.02920985221862793, 0.042002201080322266, 0.02824091911315918, 0.028715848922729492, 0.02948164939880371, 0.028994321823120117, 0.030215740203857422, 0.029191255569458008, 0.029179811477661133, 0.028284311294555664, 0.028439760208129883, 0.031086444854736328, 0.028697490692138672, 0.029597997665405273, 0.029361963272094727, 0.028464317321777344, 0.031247854232788086, 0.029204130172729492, 0.03040599822998047, 0.03568601608276367, 0.03097677230834961, 0.02924489974975586, 0.027979373931884766, 0.029233694076538086, 0.02788710594177246, 0.02768254280090332, 0.03287816047668457, 0.03187155723571777, 0.029124975204467773, 0.029724597930908203, 0.028680086135864258, 0.02904677391052246, 0.029363632202148438, 0.02868175506591797, 0.028097867965698242, 0.027370214462280273, 0.031381845474243164, 0.028300762176513672]
Trace ELOD
[55.79166666666674, 141.96226415094338, 94.5, 64.76271186440681, 102.61538461538453, 147.36936936936945, 54.05555555555554, 112.96190476190486, 140.57142857142867, 60.86956521739131, 72.64583333333326, 46.41666666666663, 58.86538461538464, 28.400000000000006, 68.65957446808511, 67.0, 115.17171717171732, 71.80769230769226, 122.16504854368941, 88.08695652173913, 140.75, 84.19999999999993, 75.39999999999998, 70.36734693877554, 116.98000000000002, 93.50847457627117, 130.47000000000003, 122.87000000000012, 81.32142857142856, 58.36363636363632, 62.9375, 87.0, 78.56363636363642, 71.60416666666663, 80.0181818181818, 62.931818181818244, 53.66666666666674, 134.28571428571433, 68.14285714285711, 128.0695652173913, 88.81632653061217, 61.415094339622556, 74.0769230769231, 128.23404255319156, 128.77777777777783, 129.1649484536083, 138.66666666666663, 77.0, 125.30303030303025, 81.8039215686274, 124.834862385321, 70.90909090909088, 69.44444444444446, 103.39999999999998, 66.0, 76.56000000000006, 65.875, 138.41121495327116, 58.43478260869563, 71.625, 44.18181818181819, 26.181818181818187, 77.9230769230769, 61.66666666666663, 36.46153846153847, 85.05555555555554, 109.41666666666674, 78.09090909090912, 65.57142857142856, 76.1875, 113.2790697674418, 68.77777777777771, 72.94117647058829, 78.21052631578948, 81.29411764705878, 138.4905660377358, 34.0, 85.03846153846155, 64.695652173913, 123.66666666666674, 111.32653061224494, 69.52307692307693, 67.66666666666663, 85.0, 79.85185185185185, 64.06666666666672, 62.67999999999995, 80.96153846153845, 71.10344827586209, 66.95121951219517, 127.95876288659792, 37.380952380952294, 81.80000000000001, 72.76086956521749, 50.07894736842104, 77.5, 24.5, 44.636363636363626, 119.8632478632478, 37.15789473684208]
PCST ELOD
[15.75, 89.41509433962267, 23.774193548387096, 20.525423728813564, 39.230769230769226, 90.25225225225222, 25.38888888888889, 69.6190476190477, 89.85714285714289, 33.78260869565217, 29.1875, 33.79166666666666, 21.61538461538464, 20.2, 18.297872340425528, 29.5, 57.44444444444457, 30.769230769230774, 69.04854368932047, 29.17391304347825, 70.83333333333348, 26.0, 25.200000000000045, 41.48979591836735, 63.710000000000036, 27.37288135593218, 97.15000000000009, 108.17000000000007, 32.178571428571445, 18.545454545454533, 32.58333333333337, 24.083333333333258, 38.19999999999999, 29.125, 17.74545454545455, 26.204545454545496, 17.0, 67.85714285714289, 22.714285714285715, 63.47826086956525, 34.40816326530609, 23.622641509433947, 26.28846153846154, 69.03191489361711, 83.87962962962956, 101.56701030927843, 66.66666666666663, 33.0, 40.24242424242425, 35.372549019607845, 72.38532110091728, 32.0, 33.18518518518516, 23.461538461538453, 35.25, 29.460000000000008, 28.274999999999977, 70.12149532710282, 24.739130434782624, 29.8125, 22.454545454545467, 0.0, 41.923076923076906, 20.333333333333258, 24.615384615384617, 19.31481481481481, 40.5, 49.18181818181819, 26.5, 24.125, 63.883720930232585, 20.44444444444443, 30.29411764705884, 21.57894736842104, 34.078431372549005, 98.46226415094338, 25.0, 34.53846153846155, 22.86956521739131, 78.33333333333326, 65.57142857142856, 36.569230769230785, 37.0, 19.285714285714278, 46.111111111111086, 34.0, 22.359999999999985, 32.88461538461539, 22.17241379310346, 28.073170731707364, 68.75257731958766, 16.904761904761898, 23.04000000000002, 24.065217391304373, 29.526315789473642, 25.388888888888914, 20.049999999999983, 19.0, 60.64957264957275, 21.210526315789465]

"""






