import ROOT as r

def fit1(entries=1000, save=False):
    # entries is the number of random samples filled into the histogram

    randomHist1 = r.TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100)
    randomHist1.Sumw2()
    generator=r.TRandom2(0)  # parameter == seed, 0->use clock

    for i in range(entries):
        randomHist1.Fill(generator.Gaus(50,10)) # params: mean, sigma

    #simple fits may be performed automatically
    r.gStyle.SetOptFit(1111) # show reduced chi2, probability, and params
    randomHist1.Fit("gaus")
    randomHist1.DrawCopy("e")  # "e" shows bin errors
    # Using DrawCopy vs Draw allows us to delete the original histogram
    # without removing it from the display.  If we save the histogran to a
    # file and close the file, it will be deleted from memory.

    # Above we used a built in function, gaus, in the fit
    # This function will be associated with the histogram
    # and may be retrieved to get parameter information
    # Refer to http://root.cern.ch/root/html/TF1.html
    # for a complete list of TF1 methods

    fitfunc = randomHist1.GetFunction("gaus")
    print("\nFit Params and errors")
    for i in range(3):
        print(f'{fitfunc.GetParameter(i):.2f} +- {fitfunc.GetParError(i):.2f}')

    print(f'Fit Probability: {fitfunc.GetProb():.2f}') # returns chi^2 p-value

    if save:
        tf.Write()
        tf.Close()

    return randomHist1
# **************************************

def fit1a(entries=1000, experiments=1000):
    randomHist1 = r.TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100)
    randomHist1.Sumw2()
    generator=r.TRandom2(0)  # parameter == seed, 0->use clock

    chi_hist = r.TH1F("Chi", "Reduced Chi2;x;frequency", 100, 0, 3)
    prob_hist = r.TH1F("Prob", "Probability;x;frequency", 100, 0, 1)
    param_hist = r.TH1F("Param", "Mean;x;frequency", 100, 40, 60)
    err_hist = r.TH1F("Err", "Error;x;frequency", 100, 0.1, 0.5)
    
    for trial in range(experiments):
        randomHist1.Reset()
        for i in range(entries):
            randomHist1.Fill(generator.Gaus(50,10)) # params: mean, sigma

    #simple fits may be performed automatically
        #r.gStyle.SetOptFit(1111) # show reduced chi2, probability, and params
        gaus=r.TF1("mygaus","[0]*exp(-(x-[1])*(x-[1])/[2]/[2])",0,100)
        gaus.SetParameter(0,randomHist1.GetMaximum())
        gaus.SetParameter(1,randomHist1.GetMean()) 
        gaus.SetParameter(2,randomHist1.GetStdDev()) 
        randomHist1.Fit("mygaus", "Q N")
        #if not result.IsValid():
         #   continue
        chi2_fit=gaus.GetChisquare()
        ndf = gaus.GetNDF()
        prob = gaus.GetProb()
        mean = gaus.GetParameter(1)
        err = gaus.GetParError(1)

        

        chi_hist.Fill(chi2_fit/ndf)
        prob_hist.Fill(prob)
        param_hist.Fill(mean)
        err_hist.Fill(err)
        

    tc = r.TCanvas("canvas", "Chi2 Graphs", 800, 600)
    tc.Divide(2,2)
    tc.cd(1)
    chi_hist.Draw()
    tc.cd(2)
    param_hist.Draw()
    tc.cd(3)
    prob_hist.Draw()
    tc.cd(4)
    err_hist.Draw()
    tc.Update()
    tc.SaveAs("result1.pdf")

    



    

if __name__ == "__main__":
    #fit1()
    fit1a()
    input("hit Enter to exit")


# example of plotting a similar histogram with error bars using numpy/matplotlib
# then using lmfit to perform the fit

#import numpy as np
#from matplotlib import pyplot as plt
#
#entries=1000
#
#vals=np.random.normal(loc=50, scale=10, size=entries)
#y,binEdges=np.histogram(vals, bins=50, range=(1,100))
#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#ey         = np.sqrt(y)
#width      = binEdges[1]-binEdges[0]
#plt.bar(bincenters, y, width=width, color='r', yerr=ey)
#plt.show(block=False)
