import ROOT
from math import exp

def createCanvasPads():
    c = ROOT.TCanvas("c", "canvas", 800, 700)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.30, 1, 1.0)
    pad1.SetBottomMargin(0.02)  # joins upper and lower plot
    pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.29)
    pad2.SetTopMargin(0.01)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()
 
    return c, pad1, pad2

#Working muon file
#filename = "/work/jheitkoetter/Z_early_Run3/input_histos/earlyRun3_crown_2022_mm_run.root"


#Working file
pathname = "/work/dseyler/WpTAnalysis_fromJulius/OUTPUTS_allV10Runs_prelim"
filename = pathname + "/root/QCD/qcdshape_extrapolated_muminus_met_WpT_bin0_muminus.root"

file = ROOT.TFile(filename)

all_hist_names = [key.GetName() for key in file.GetListOfKeys()]
hist_names = []
for name in all_hist_names:
    if "QCD" not in name and "par0" not in name:
        hist_names.append(name)

#hist = file.Get("data#mm-355680-scale_1_pb#Nominal#m_vis")
for i in hist_names:
    hist_name = i
    c1, pad1, pad2 = createCanvasPads()
    """
    c1 = ROOT.TCanvas( 'c1', 'The Fit Canvas', 200, 10, 700, 500 )
    c1.SetGridx()
    c1.SetGridy()
    c1.GetFrame().SetFillColor( 21 )
    c1.GetFrame().SetBorderMode(-1 )
    c1.GetFrame().SetBorderSize( 5 )
    c1.SetTitle("Met QCD Background Pol1")
    """

    hist = file.Get(hist_name)
    hist.Print()
    print(hist.Integral())

    x = ROOT.RooRealVar("x", "x", 0, 120)
    arglist = ROOT.RooArgList(x)
    data = ROOT.RooDataHist("data", "data", arglist, hist)

    #Setting up variable fit parameters
    sigma_0 = ROOT.RooRealVar("sigma_0", "sigma_0", 25, 0, 1000)
    sigma_1 = ROOT.RooRealVar("sigma_1", "sigma_1", 0, 0, 500)
    sigma_2 = ROOT.RooRealVar("sigma_2", "sigma_2", 0, 0, 500)

    #Adding pdf
    genpdf = ROOT.RooGenericPdf("genpdf", "genpdf", "x*exp(-x**2/(sigma_0 + sigma_1*x**2 + sigma_2*x))", ROOT.RooArgList(x, sigma_0, sigma_1, sigma_2))

    #Fitting FunctionZfitResult = conv_pdf.fitTo(data)
    pad1.cd()
    fitline = genpdf.fitTo(data)
    frame = x.frame()
    data.plotOn(frame)
    genpdf.plotOn(frame)

    frame.GetYaxis().SetLabelSize(0.04)


    frame.Draw()
    """
    frame.GetYaxis().SetLabelSize(0.0)
    axis = ROOT.TGaxis(0, 100000, 0, 770000, 100000, 770000, 510, "")
    axis.SetLabelFont(43)
    axis.SetLabelSize(20)
    axis.Draw()
    """
    #Axes
    lepton_types = ["muplus", "muminus", "eplus", "eminus"]
    lepton = ""
    for i in lepton_types:
        if i in hist_name and lepton == "":
            lepton = i
    fit_types = ["pol1", "pol2"]
    fit_type = ""
    for i in fit_types:
        if i in hist_name and fit_type == "":
            fit_type = i
    modifier = ""
    if "prelim" in filename:
        modifier = "prelim iso mean"
    elif "trueMean" in filename:
        modifier = "true iso mean"
    elif "eta_binned" in filename:
        if "bin2" in hist_name:
            modifier = "endcap"
        elif "bin1" in hist_name:
            modifier = "barrel"
        else:
            modifier = "all eta"


    Title = "WpT QCD background " + lepton + " " + fit_type + " " + modifier
    frame.SetTitle(Title)
    frame.GetYaxis().SetTitle("Count")
    frame.GetXaxis().SetLabelSize(0.0)




    #Text

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.DrawText(0.5,0.7, "sigma_0 = %.4f "%(sigma_0.getVal()) + " +/- %.4f"%(sigma_0.getError()))
    latex.DrawText(0.5,0.65, "sigma_1 = 0%.4f "%(sigma_1.getVal()) + " +/- %.4f"%(sigma_1.getError()))
    latex.DrawText(0.5,0.6, "sigma_2 = %.4f "%(sigma_2.getVal()) + " +/- %.4f"%(sigma_2.getError()))
    latex.DrawText(0.5, 0.55, "Chi2/ndof = %.4f "%(21*frame.chiSquare())+"/21 " + "= %.4f"%(frame.chiSquare()))

    # -------------------------------------------------------
 
    # Construct a histogram with the residuals of the data w.r.t. the curve
    hresid = frame.residHist()
 
    # Construct a histogram with the pulls of the data w.r.t the curve
    hpull = frame.pullHist()
 
    # Create a frame to draw the residual distribution and add the
    # distribution to the frame
    frame2 = x.frame()#Title="Residual Distribution")
    frame2.addPlotable(hresid, "P")
    
    # Create a frame to draw the pull distribution and add the distribution to
    # the frame
    frame3 = x.frame()#Title="Pull Distribution")
    frame3.addPlotable(hpull, "P")


    frame3.GetXaxis().SetTitle("Met [GeV]")
    #frame3.GetYaxis().SetTitle("Pull")
    frame3.GetYaxis().SetLabelSize(0.12)
    frame3.GetXaxis().SetLabelSize(0.12)

    frame3.SetTitle("")
    frame3.GetXaxis().SetTitleSize(0.1)
    

    pad2.cd()
    frame3.Draw()
    #print(frac.getVal())

    print(hist_name)
    print("sigma_0: ", sigma_0.getVal(), " +- ", sigma_0.getError())
    print("sigma_1: ", sigma_1.getVal(), " +- ", sigma_1.getError())
    print("sigma_2: ", sigma_2.getVal(), " +- ", sigma_2.getError())
    print("Chi2:", frame.chiSquare())
    print(hist_names)

    #Legend
    '''
    legend = ROOT.TLegend(0.7,0.6,0.85,0.75)
    legend.AddEntry(hist, "Data")
    legend.AddEntry(fitline, "Fit")
    legend.SetLineWidth(2)
    legend.Draw("same")
    '''

   #Descriptive name for the output file
    modify_name = ""
    if "prelim" in filename:
        modify_name = "prelim"
    elif "trueMean" in filename:
        modify_name = "trueMean"
    elif "eta_binned" in filename:
        if "bin2" in hist_name:
            modify_name = "endcap"
        elif "bin1" in hist_name:
            modify_name = "barrel"
        else:
            modify_name = "allEta"


    c1.Update()
    c1.SaveAs("/work/dseyler/WpTAnalysis_fromJulius/plots/png/Fit_" + lepton + "_" + fit_type + "_" + modify_name + ".png")
    c1.SaveAs("/work/dseyler/WpTAnalysis_fromJulius/plots/pdf/Fit_" + lepton + "_" + fit_type + "_" + modify_name + ".pdf")
    c1.SaveAs("/work/dseyler/WpTAnalysis_fromJulius/plots/root/Fit_" + lepton + "_" + fit_type + "_" + modify_name + ".root")