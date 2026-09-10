// Three-way comparison of the profile likelihood in each oscillation parameter.
//
//   Profiling  outputs/profile/  one independent gundamFitter per point, the
//                                scanned parameter isFrozen, everything else
//                                refit. RootMinimizer, strategy 1, tol 1E-2.
//   1D grid    outputs/grid1d/   GridSearch on that parameter alone, the other
//                                two oscillation parameters left free to Migrad.
//                                strategy 0, tolerance 1, warm start.
//   3D grid    outputs/grid/     the full grid, profiled by taking the minimum
//                                over the two axes not being plotted.
//
//   root -l -b -q 'compare3way.C("outputs/profile","outputs/grid1d","outputs/grid","outputs/compare3way.root")'


const double kD2R = TMath::Pi()/180.;

struct Par {
    const char* key;        // file stem in outputs/profile and outputs/grid1d
    const char* axisSub;    // identifies the axis inside AxisName
    int gridIndex;          // position in the 3D grid's GridValues[]
    double scale;           // multiply stored value by this before keying
    int round;              // decimal places used to key the map
    const char* title;
    const char* unit;
};

// keys are rounded because a grid value and a scan value that mean the same
// point do not agree bit for bit
double key(double v, int nd){ double p = TMath::Power(10,nd); return TMath::Nint(v*p)/p; }

TFile* gOut{nullptr};
std::string gProfileDir, gGrid1dDir, gGridDir;

std::map<double,double> readProfile(const Par& p){
    std::map<double,double> out;
    TSystemDirectory dir("d", gProfileDir.c_str());
    auto* files = dir.GetListOfFiles();
    if(!files) return out;
    for(auto o : *files){
        TString n = ((TSystemFile*)o)->GetName();
        if(!n.EndsWith(".root")) continue;
        if(!n.Contains(Form("_With_%s_", p.key))) continue;

        TFile f(TString(gProfileDir) + "/" + n);
        auto t = f.Get<TTree>("FitterEngine/postFit/bestFitStats");
        auto s = f.Get<TNamed>("FitterEngine/postFit/parState_TNamed");
        if(!t || !s) continue;

        double llh; t->SetBranchAddress("totalLikelihoodAtBestFit",&llh); t->GetEntry(0);

        // the frozen value is not stored as a number anywhere, so take it from
        // the injected parameter state that the job was given
        TString js(s->GetTitle());
        TString tag = p.axisSub;
        int i = js.Index(tag);
        if(i < 0) continue;
        int v = js.Index("\"value\"", i);
        if(v < 0) continue;
        TString num = js(v+8, 40);
        num.ReplaceAll(":",""); num.ReplaceAll(" ","");
        double val = TString(num(0, num.First(reinterpret_cast<const char*>("\n")[0]))).Atof();
        if(val == 0) val = atof(num.Data());

        out[ key(val*p.scale, p.round) ] = llh;
    }
    return out;
}

std::map<double,double> read1D(const Par& p){
    std::map<double,double> out;
    TFile f(Form("%s/grid1d_%s.root", gGrid1dDir.c_str(), p.key));
    auto g = f.Get<TTree>("FitterEngine/fit/gridTree");
    if(!g) return out;
    double llh, gv[8];
    g->SetBranchAddress("LLH",&llh);
    g->SetBranchAddress("GridValues",gv);
    for(Long64_t i=0;i<g->GetEntries();i++){ g->GetEntry(i);
        out[ key(gv[0]*p.scale, p.round) ] = llh; }
    return out;
}

std::map<double,double> read3D(const Par& p){
    std::map<double,double> out;
    TChain c("FitterEngine/fit/gridTree");
    for(int i=0;i<51;i++) c.Add(Form("%s/grid_%d.root", gGridDir.c_str(), i));
    double llh, gv[3];
    c.LoadTree(0);
    // GridValues replaced Grided; read whichever these files have
    const char* vBranch = c.GetBranch("GridValues") ? "GridValues" : "Grided";
    c.SetBranchAddress("LLH",&llh);
    c.SetBranchAddress(vBranch,gv);
    for(Long64_t i=0;i<c.GetEntries();i++){ c.GetEntry(i);
        double k = key(gv[p.gridIndex]*p.scale, p.round);
        if(!out.count(k) || llh < out[k]) out[k] = llh; }
    return out;
}

double minOf(const std::map<double,double>& m){
    double b=1e30; for(auto& kv : m) if(kv.second<b) b=kv.second; return b;
}

void draw(const Par& p,
          std::map<double,double>& m1, std::map<double,double>& m2, std::map<double,double>& m3){
    double b1=minOf(m1), b2=minOf(m2), b3=minOf(m3);

    auto g1 = new TGraph(), g2 = new TGraph(), g3 = new TGraph();
    for(auto& kv : m1) g1->SetPoint(g1->GetN(), kv.first, kv.second-b1);
    for(auto& kv : m2) g2->SetPoint(g2->GetN(), kv.first, kv.second-b2);
    for(auto& kv : m3) g3->SetPoint(g3->GetN(), kv.first, kv.second-b3);

    auto cv = new TCanvas("cv","",900,700);
    auto pad1 = new TPad("pad1","",0,0.30,1,1);   pad1->SetBottomMargin(0.02);
    auto pad2 = new TPad("pad2","",0,0.00,1,0.30);pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.30);
    pad1->SetLeftMargin(0.12); pad2->SetLeftMargin(0.12);
    pad1->SetGridx(); pad1->SetGridy(); pad2->SetGridx(); pad2->SetGridy();
    pad1->Draw(); pad2->Draw();

    double xlo=1e30, xhi=-1e30;
    for(auto& kv : m1){ xlo=TMath::Min(xlo,kv.first); xhi=TMath::Max(xhi,kv.first); }

    // whole curve, no cap. headroom leaves the top strip clear for the legend
    double ymax=0;
    for(int i=0;i<g1->GetN();i++) ymax=TMath::Max(ymax,g1->GetY()[i]);
    for(int i=0;i<g2->GetN();i++) ymax=TMath::Max(ymax,g2->GetY()[i]);
    for(int i=0;i<g3->GetN();i++) ymax=TMath::Max(ymax,g3->GetY()[i]);
    ymax *= 1.45;

    pad1->cd();
    auto fr = pad1->DrawFrame(xlo, 0, xhi, ymax);
    fr->GetYaxis()->SetTitle("#Delta#chi^{2}");
    fr->GetYaxis()->SetTitleSize(0.055);  fr->GetYaxis()->SetTitleOffset(1.05);
    fr->GetYaxis()->SetLabelSize(0.045);
    fr->GetXaxis()->SetLabelSize(0);

    g1->SetLineColor(kBlack);    g1->SetLineStyle(1); g1->SetLineWidth(4); g1->Draw("L SAME");
    g2->SetLineColor(kRed+1);    g2->SetLineStyle(2); g2->SetLineWidth(3); g2->Draw("L SAME");
    g3->SetLineColor(kAzure+2);  g3->SetLineStyle(3); g3->SetLineWidth(3); g3->Draw("L SAME");

    auto lg = new TLegend(0.60,0.70,0.90,0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
    lg->AddEntry(g1,"Profiling","L");
    lg->AddEntry(g2,"1D grid","L");
    lg->AddEntry(g3,"3D grid","L");
    lg->Draw();

    auto tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.040);
    tx->DrawLatex(0.15,0.90, Form("%s", p.title));
    tx->SetTextSize(0.033);
    tx->DrawLatex(0.15,0.845, Form("minima  %.4f / %.4f / %.4f", b1, b2, b3));

    pad2->cd();
    auto fr2 = pad2->DrawFrame(xlo, -0.5, xhi, 0.5);
    fr2->GetXaxis()->SetTitle(Form("%s%s", p.title, p.unit));
    fr2->GetYaxis()->SetTitle("#minus Profiling");
    fr2->GetYaxis()->SetTitleSize(0.10); fr2->GetYaxis()->SetTitleOffset(0.45);
    fr2->GetYaxis()->SetLabelSize(0.08); fr2->GetYaxis()->SetNdivisions(505);
    fr2->GetXaxis()->SetTitleSize(0.10); fr2->GetXaxis()->SetLabelSize(0.09);

    auto d2 = new TGraph(), d3 = new TGraph();
    for(auto& kv : m1){
        double a = kv.second-b1;
        if(m2.count(kv.first)) d2->SetPoint(d2->GetN(), kv.first, (m2[kv.first]-b2)-a);
        if(m3.count(kv.first)) d3->SetPoint(d3->GetN(), kv.first, (m3[kv.first]-b3)-a);
    }
    d2->SetLineColor(kRed+1);   d2->SetLineStyle(2); d2->SetLineWidth(3); d2->Draw("L SAME");
    d3->SetLineColor(kAzure+2); d3->SetLineStyle(3); d3->SetLineWidth(3); d3->Draw("L SAME");

    cv->SetName(Form("compare3way_%s", p.key));
    cv->SetTitle(p.title);
    gOut->cd();
    cv->Write();
}

void report(const Par& p,
            std::map<double,double>& m1, std::map<double,double>& m2, std::map<double,double>& m3){
    double b1=minOf(m1), b2=minOf(m2), b3=minOf(m3);
    printf("\n=== %s ===\n", p.title);
    printf("points: Profiling %zu   1D grid %zu   3D grid %zu\n", m1.size(), m2.size(), m3.size());
    printf("absolute minimum: Profiling %.5f   1D grid %.5f   3D grid %.5f\n", b1, b2, b3);

    // the absolute deviation scales with how high the curve goes, so quote the
    // span too: a 1.4 discrepancy means something different on a curve of 160
    // than on one of 3000
    double span=0; for(auto& kv : m1) span = TMath::Max(span, kv.second-b1);
    printf("curve spans 0 to %.1f in dchi2\n", span);

    double s2=0,s3=0,x2=0,x3=0,at2=0,at3=0; int n2=0,n3=0;
    for(auto& kv : m1){
        double a = kv.second-b1;
        if(m2.count(kv.first)){ double d=(m2[kv.first]-b2)-a; s2+=fabs(d); if(fabs(d)>x2){x2=fabs(d);at2=kv.first;} n2++; }
        if(m3.count(kv.first)){ double d=(m3[kv.first]-b3)-a; s3+=fabs(d); if(fabs(d)>x3){x3=fabs(d);at3=kv.first;} n3++; }
    }
    if(n2) printf("1D grid vs Profiling: %d shared points, mean |d| %.4f, max %.4f at %.2f  (%.3f%% of span)\n",
                  n2, s2/n2, x2, at2, 100*x2/span);
    else   printf("1D grid vs Profiling: no shared points\n");
    if(n3) printf("3D grid vs Profiling: %d shared points, mean |d| %.4f, max %.4f at %.2f  (%.3f%% of span)\n",
                  n3, s3/n3, x3, at3, 100*x3/span);
    else   printf("3D grid vs Profiling: no shared points\n");
}

void compare3way(std::string profile_dir, std::string grid1d_dir, std::string grid_dir, std::string out_file){
    gProfileDir = profile_dir; gGrid1dDir = grid1d_dir; gGridDir = grid_dir;
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetLegendFont(42);
    TGaxis::SetMaxDigits(3);
    // frames made by DrawFrame would otherwise be written into gOut alongside
    // the canvases
    TH1::AddDirectory(kFALSE);

    gOut = new TFile(out_file.c_str(), "RECREATE");

    std::vector<Par> pars = {
        { "dcp",   "PMNS_DELTA_CP",              1, 1./kD2R, 1, "#delta_{CP}",           " [deg]" },
        { "sin23", "PMNS_SIN_SQUARED_23",        0, 1.,      4, "sin^{2}#theta_{23}",    ""       },
        { "dm32",  "PMNS_DELTA_MASS_SQUARED_32", 2, 1e5,     2, "#Deltam^{2}_{32}",      " [10^{-5} eV^{2}]" },
    };

    for(auto& p : pars){
        auto m1 = readProfile(p);
        auto m2 = read1D(p);
        auto m3 = read3D(p);
        report(p, m1, m2, m3);
        if(!m1.empty() && !m2.empty() && !m3.empty()) draw(p, m1, m2, m3);
    }

    gOut->Close();
    printf("\nwrote %s\n", out_file.c_str());
}
