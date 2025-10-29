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


def fit1b(entries=1000, experiments=1000):
    randomHist1 = r.TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100)
    randomHist1.Sumw2()
    generator=r.TRandom2(0)  # parameter == seed, 0->use clock

    chi_hist = r.TH1F("Param", "Mean Chi2;x;frequency", 100, 45, 55)
    nll_hist = r.TH1F("Err", "Mean NLL;x;frequency", 100, 45, 55)
    
    for trial in range(experiments):
        randomHist1.Reset()
        for i in range(entries):
            randomHist1.Fill(generator.Gaus(50,10)) # params: mean, sigma

    #simple fits may be performed automatically
        #r.gStyle.SetOptFit(1111) # show reduced chi2, probability, and params
        gaus_chi=r.TF1("gaus_chi","[0]*exp(-(x-[1])*(x-[1])/[2]/[2])",0,100)
        gaus_chi.SetParameter(0,randomHist1.GetMaximum())
        gaus_chi.SetParameter(1,randomHist1.GetMean()) 
        gaus_chi.SetParameter(2,randomHist1.GetStdDev()) 
        randomHist1.Fit("gaus_chi", "Q N")

        gaus_ll=r.TF1("gaus_ll","[0]*exp(-(x-[1])*(x-[1])/[2]/[2])",0,100)
        gaus_ll.SetParameter(0,randomHist1.GetMaximum())
        gaus_ll.SetParameter(1,randomHist1.GetMean()) 
        gaus_ll.SetParameter(2,randomHist1.GetStdDev()) 
        randomHist1.Fit("gaus_ll", "Q N L")

        chi = gaus_chi.GetParameter(1)
        ll = gaus_ll.GetParameter(1)
        

        

        chi_hist.Fill(chi)
        nll_hist.Fill(ll)
        

    tc = r.TCanvas("canvas", "Chi2 Graphs", 800, 600)
    tc.Divide(2,1)
    tc.cd(1)
    chi_hist.Draw()
    tc.cd(2)
    nll_hist.Draw()
    tc.Update()
    tc.SaveAs("result2.pdf")


if __name__ == "__main__":
    #fit1()
    fit1b()
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
