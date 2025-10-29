# run with: python readhist.py

import ROOT as r
import numpy as np
import matplotlib.pyplot as plt

def NegLogLik(hist, fit_func):
    count  = 0
    for i in range(1, hist.GetNbinsX()+1):
        xi = hist.GetBinCenter(i)
        lamda = fit_func.Eval(xi)
        yi = hist.GetBinContent(i)
        count += yi*np.log(lamda)-lamda-r.TMath.LnGamma(yi+1)
    return -2*count
        

file25 = r.TFile("histo25.root")
h25 = file25.Get("randomHist1")

generator=r.TRandom2(0)
nll_hist = r.TH1F("nll", "nll", 100, 50, 150)

gaus=r.TF1("gaus","[0]*exp(-(x-[1])*(x-[1])/[2]/[2])",0,100)
gaus.SetParameter(0,h25.GetMaximum())
gaus.SetParameter(1,h25.GetMean()) 
gaus.SetParameter(2,h25.GetStdDev()) 
h25.Fit("gaus", "Q N L")

amplitude = gaus.GetParameter(0)
mean = gaus.GetParameter(1)
sigma = gaus.GetParameter(2)

nll_real = NegLogLik(h25,gaus)

print(f"The NLL value is {nll_real:.3f}")

exp_hist = r.TH1F("exp", "exp", h25.GetNbinsX(), 0, 100)
gaus_exp=r.TF1("gaus_exp","[0]*exp(-(x-[1])*(x-[1])/[2]/[2])",0,100)

for i in range(1000):
    exp_hist.Reset()

    for j in range(25):
        exp_hist.Fill(generator.Gaus(50,10))

    gaus_exp.SetParameter(0,exp_hist.GetMaximum())
    gaus_exp.SetParameter(1,exp_hist.GetMean()) 
    gaus_exp.SetParameter(2,exp_hist.GetStdDev()) 
    exp_hist.Fit("gaus_exp", "Q N L")

    nll_exp = NegLogLik(exp_hist, gaus_exp)
    nll_hist.Fill(nll_exp)


tc = r.TCanvas("NLL", "NLL")
nll_hist.Draw()
line = r.TLine(nll_real, 0, nll_real, nll_hist.GetMaximum())
line.SetLineColor(r.kRed)
line.SetLineWidth(2)
line.Draw("same")
tc.Update()
tc.SaveAs("result3.pdf")



count = 0
for i in range(1, nll_hist.GetNbinsX()+1):
    if nll_hist.GetBinCenter(i) + 1 <= nll_real:
        count += nll_hist.GetBinContent(i)

print(f"The P value is {count/1000:.3f}")


means = np.linspace(-4.25,4.25)+mean
nlls = []
gaus_cont=r.TF1("gaus_cont","[0]*exp(-(x-[1])*(x-[1])/[2]/[2])",0,100)
gaus_cont.SetParameter(0, amplitude)
gaus_cont.SetParameter(2, sigma)
for i in means:
    gaus_cont.SetParameter(1, i)
    nlls.append(NegLogLik(h25, gaus_cont))

fig = plt.figure()
plt.plot(means, nlls)
plt.xlabel("Mean Values")
plt.ylabel("NLL Values")
plt.grid()
plt.show()




    

input("hit return to exit")

