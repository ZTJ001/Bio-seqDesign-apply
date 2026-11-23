#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""目标蛋白 突变效应预测模型的特征构造器（公开版）。

仅含标准 BLOSUM62 与氨基酸物化常数，不含任何实验数据或数据来源信息。
53 维特征顺序必须与模型训练时完全一致，不可调整。
"""
import math

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
CHARGE = {**{a: 1 for a in "KRH"}, **{a: -1 for a in "DE"}}
KD = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
VOL = {"A":88.6,"C":108.5,"D":111.1,"E":138.4,"F":189.9,"G":60.1,"H":153.2,"I":166.7,"K":168.6,"L":166.7,"M":162.9,"N":114.1,"P":112.7,"Q":143.8,"R":173.4,"S":89.0,"T":116.1,"V":140.0,"W":227.8,"Y":193.6}

BLOSUM62 = {
"AAA":4,"AAC":0,"AAD":-2,"AAE":-1,"AAF":-2,"AAG":0,"AAH":-2,"AAI":-1,"AAK":-1,"AAL":-1,"AAM":-1,"AAN":-2,"AAP":-1,"AAQ":-1,"AAR":-2,"AAS":1,"AAT":0,"AAV":0,"AAW":-3,"AAY":-2,
"CCA":0,"CCC":9,"CCD":-3,"CCE":-4,"CCF":-2,"CCG":-3,"CCH":-3,"CCI":-1,"CCK":-3,"CCL":-1,"CCM":-1,"CCN":-3,"CCP":-3,"CCQ":-3,"CCR":-3,"CCS":-1,"CCT":-1,"CCV":-1,"CCW":-2,"CCY":-2,
"DDA":-2,"DDC":-3,"DDD":6,"DDE":2,"DDF":-3,"DDG":-1,"DDH":-1,"DDI":-3,"DDK":-1,"DDL":-4,"DDM":-3,"DDN":1,"DDP":-1,"DDQ":0,"DDR":-2,"DDS":0,"DDT":-1,"DDV":-3,"DDW":-4,"DDY":-3,
"EEA":-1,"EEC":-4,"EED":2,"EEE":5,"EEF":-3,"EEG":-2,"EEH":-2,"EEI":-3,"EEK":1,"EEL":-3,"EEM":-2,"EEN":0,"EEP":-1,"EEQ":2,"EER":0,"EES":0,"EET":-1,"EEV":-2,"EEW":-3,"EEY":-2,
"FFA":-2,"FFC":-2,"FFD":-3,"FFE":-3,"FFF":6,"FFG":-2,"FFH":-1,"FFI":0,"FFK":-3,"FFL":1,"FFM":0,"FFN":-2,"FFP":-4,"FFQ":-3,"FFR":-2,"FFS":-2,"FFT":-2,"FFV":0,"FFW":1,"FFY":3,
"GGA":0,"GGC":-3,"GGD":-1,"GGE":-2,"GGF":-2,"GGG":6,"GGH":-2,"GGI":-4,"GGK":-2,"GGL":-4,"GGM":-3,"GGN":0,"GGP":-2,"GGQ":-2,"GGR":-2,"GGS":0,"GGT":-1,"GGV":-3,"GGW":-2,"GGY":-3,
"HHA":-2,"HHC":-3,"HHD":-1,"HHE":-2,"HHF":-1,"HHG":-2,"HHH":8,"HHI":-3,"HHK":-2,"HHL":-3,"HHM":-2,"HHN":1,"HHP":-2,"HHQ":2,"HHR":0,"HHS":-1,"HHT":-2,"HHV":-3,"HHW":-2,"HHY":2,
"IIA":-1,"IIC":-1,"IID":-3,"IIE":-3,"IIF":0,"IIG":-4,"IIH":-3,"III":4,"IIK":-3,"IIL":2,"IIM":1,"IIN":-3,"IIP":-3,"IIQ":-3,"IIR":-3,"IIS":-2,"IIT":-1,"IIV":3,"IIW":-3,"IIY":-1,
"KKA":-1,"KKC":-3,"KKD":-1,"KKE":1,"KKF":-3,"KKG":-2,"KKH":-2,"KKI":-3,"KKK":5,"KKL":-2,"KKM":-1,"KKN":0,"KKP":-1,"KKQ":2,"KKR":2,"KKS":0,"KKT":-1,"KKV":-2,"KKW":-3,"KKY":-2,
"LLA":-1,"LLC":-1,"LLD":-4,"LLE":-3,"LLF":1,"LLG":-4,"LLH":-3,"LLI":2,"LLK":-2,"LLL":4,"LLM":2,"LLN":-3,"LLP":-3,"LLQ":-2,"LLR":-2,"LLS":-2,"LLT":-1,"LLV":1,"LLW":-2,"LLY":-1,
"MMA":-1,"MMC":-1,"MMD":-3,"MME":-2,"MMF":0,"MMG":-3,"MMH":-2,"MMI":1,"MMK":-1,"MML":2,"MMM":5,"MMN":-2,"MMP":-2,"MMQ":0,"MMR":-1,"MMS":-1,"MMT":-1,"MMV":1,"MMW":-1,"MMY":-1,
"NNA":-2,"NNC":-3,"NND":1,"NNE":0,"NNF":-2,"NNG":0,"NNH":1,"NNI":-3,"NNK":0,"NNL":-3,"NNM":-2,"NNN":6,"NNP":-2,"NNQ":0,"NNR":0,"NNS":1,"NNT":0,"NNV":-3,"NNW":-4,"NNY":-2,
"PPA":-1,"PPC":-3,"PPD":-1,"PPE":-1,"PPF":-4,"PPG":-2,"PPH":-2,"PPI":-3,"PPK":-1,"PPL":-3,"PPM":-2,"PPN":-2,"PPP":7,"PPQ":-1,"PPR":-2,"PPS":-1,"PPT":-1,"PPV":-2,"PPW":-4,"PPY":-3,
"QQA":-1,"QQC":-3,"QQD":0,"QQE":2,"QQF":-3,"QQG":-2,"QQH":2,"QQI":-3,"QQK":2,"QQL":-2,"QQM":0,"QQN":0,"QQP":-1,"QQQ":5,"QQR":1,"QQS":-1,"QQT":-1,"QQV":-2,"QQW":-2,"QQY":-1,
"RRA":-2,"RRC":-3,"RRD":-2,"RRE":0,"RRF":-2,"RRG":-2,"RRH":0,"RRI":-3,"RRK":2,"RRL":-2,"RRM":-1,"RRN":0,"RRP":-2,"RRQ":1,"RRR":5,"RRS":-1,"RRT":-1,"RRV":-3,"RRW":-3,"RRY":-1,
"SSA":1,"SSC":-1,"SSD":0,"SSE":0,"SSF":-2,"SSG":0,"SSH":-1,"SSI":-2,"SSK":0,"SSL":-2,"SSM":-1,"SSN":1,"SSP":-1,"SSQ":-1,"SSR":-1,"SSS":4,"SST":1,"SSV":-2,"SSW":-3,"SSY":-2,
"TTA":0,"TTC":-1,"TTD":-1,"TTE":-1,"TTF":-2,"TTG":-1,"TTH":-2,"TTI":-1,"TTK":-1,"TTL":-1,"TTM":-1,"TTN":0,"TTP":-1,"TTQ":-1,"TTR":-1,"TTS":1,"TTT":5,"TTV":0,"TTW":-2,"TTY":-2,
"VVA":0,"VVC":-1,"VVD":-3,"VVE":-2,"VVF":0,"VVG":-3,"VVH":-3,"VVI":3,"VVK":-2,"VVL":1,"VVM":1,"VVN":-3,"VVP":-2,"VVQ":-2,"VVR":-3,"VVS":-2,"VVT":0,"VVV":4,"VVW":-3,"VVY":-1,
"WWA":-3,"WWC":-2,"WWD":-4,"WWE":-3,"WWF":1,"WWG":-2,"WWH":-2,"WWI":-3,"WWK":-3,"WWL":-2,"WWM":-1,"WWN":-4,"WWP":-4,"WWQ":-2,"WWR":-3,"WWS":-3,"WWT":-2,"WWV":-3,"WWW":11,"WWY":2,
"YYA":-2,"YYC":-2,"YYD":-3,"YYE":-2,"YYF":3,"YYG":-3,"YYH":2,"YYI":-1,"YYK":-2,"YYL":-1,"YYM":-1,"YYN":-2,"YYP":-3,"YYQ":-1,"YYR":-1,"YYS":-2,"YYT":-2,"YYV":-1,"YYW":2,"YYY":7,
}

CLASS_OF = lambda a: ("正电" if a in "KRH" else "负电" if a in "DE" else
                      "疏水" if a in "AVILMFWYC" else "极性" if a in "STNQ" else "特殊")


def featurize(wt, site, mut):
    """按训练时的顺序构造 53 维特征向量。"""
    v = [site/133.0, math.sin(site/133*2*math.pi), math.cos(site/133*2*math.pi)]
    ow = [0.0]*20; ow[AA_IDX[wt]] = 1.0
    om = [0.0]*20; om[AA_IDX[mut]] = 1.0
    v += ow + om
    v += [BLOSUM62.get(wt+wt+mut, 0)/4.0,
          1.0 if CLASS_OF(wt) == CLASS_OF(mut) else 0.0,
          (CHARGE.get(mut,0)-CHARGE.get(wt,0))/2.0,
          (VOL[mut]-VOL[wt])/200.0,
          (KD[mut]-KD[wt])/10.0,
          1.0 if mut == "P" else 0.0, 1.0 if mut == "G" else 0.0, 1.0 if mut == "W" else 0.0,
          VOL[wt]/250.0, KD[wt]/5.0]
    return v
