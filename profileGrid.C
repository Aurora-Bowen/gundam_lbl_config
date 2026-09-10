// Profile a GridSearch scan of any dimensionality: every axis on its own, and
// every pair of axes. Axes, node values and parameter names all come from
// axisTree, so nothing here is specific to a particular grid.
//
//   root -l -b -q 'profileGrid.C("outputs/grid/grid_*.root","outputs/grid_profiles.root")'
//
// Contours are drawn at raw Delta(chi2) values. Turning those into coverage
// needs Feldman-Cousins, so no sigma or percentage labels appear anywhere.


const double kContours[] = { 2.30, 4.61, 11.83 };
const int    kNContours  = 3;

// dCP is stored in radians but is unreadable that way
const double kRad2Deg = 180./TMath::Pi();

// axis values are rescaled rather than left to TGaxis::SetMaxDigits, which is a
// global read at paint time and so does not survive into a saved canvas
double axisScale(const std::string& n){
    if( n=="PMNS_DELTA_CP" )              return kRad2Deg;
    if( n=="PMNS_DELTA_MASS_SQUARED_32" ) return 1e3;
    if( n=="PMNS_DELTA_MASS_SQUARED_21" ) return 1e5;
    return 1.;
}

struct Axis {
    std::string name;            // full title from axisTree
    std::string label;           // for the plot
    std::vector<double> values;  // node values, already scaled
    std::vector<double> edges;   // bin edges built from the node values
    double scale{1.};
};

// ROOT needs bin edges, and the node values are bin centres. Midpoints between
// neighbours, with the outer two extrapolated. Handles a non-uniform axis.
std::vector<double> edgesFrom(const std::vector<double>& v){
    std::vector<double> e;
    if( v.size() == 1 ){ e = { v[0]-0.5, v[0]+0.5 }; return e; }
    e.push_back( v[0] - 0.5*(v[1]-v[0]) );
    for( size_t i=1;i<v.size();i++ ) e.push_back( 0.5*(v[i-1]+v[i]) );
    e.push_back( v.back() + 0.5*(v.back()-v[v.size()-2]) );
    return e;
}

// the parameter's short name, dropping the set and the "#N_" ordering prefix
std::string shortName(const std::string& full){
    size_t s = full.rfind('/');
    std::string n = (s==std::string::npos) ? full : full.substr(s+1);
    size_t u = n.find('_');
    if( u!=std::string::npos && n[0]=='#' ) n = n.substr(u+1);
    return n;
}

std::string prettyLabel(const std::string& n){
    if( n=="PMNS_DELTA_CP" )              return "#delta_{CP} [deg]";
    if( n=="PMNS_SIN_SQUARED_23" )        return "sin^{2}#theta_{23}";
    if( n=="PMNS_SIN_SQUARED_13" )        return "sin^{2}#theta_{13}";
    if( n=="PMNS_SIN_SQUARED_12" )        return "sin^{2}#theta_{12}";
    if( n=="PMNS_DELTA_MASS_SQUARED_32" ) return "#Deltam^{2}_{32} [10^{-3} eV^{2}]";
    if( n=="PMNS_DELTA_MASS_SQUARED_21" ) return "#Deltam^{2}_{21} [10^{-5} eV^{2}]";
    return n;
}

void profileGrid(std::string filePattern, std::string outPath){

    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetLegendFont(42);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);
    TGaxis::SetMaxDigits(3);
    TH1::AddDirectory(kFALSE);

    // ---- axis definitions -------------------------------------------------
    TChain ax("FitterEngine/fit/axisTree");
    if( ax.Add(filePattern.c_str()) == 0 ){ printf("no files match %s\n", filePattern.c_str()); return; }

    std::vector<std::string>* axName{nullptr};
    std::vector<std::vector<double>>* axValues{nullptr};
    ax.SetBranchAddress("AxisName", &axName);
    ax.SetBranchAddress("AxisValues", &axValues);
    ax.GetEntry(0);

    std::vector<Axis> axes;
    for( size_t i=0;i<axName->size();i++ ){
        Axis a;
        a.name  = shortName( axName->at(i) );
        a.scale = axisScale( a.name );
        a.label = prettyLabel( a.name );
        for( double v : axValues->at(i) ) a.values.push_back( v*a.scale );
        std::sort( a.values.begin(), a.values.end() );
        a.edges = edgesFrom( a.values );
        axes.push_back(a);
    }
    const int nAx = int(axes.size());
    printf("%d axes over %d files\n", nAx, ax.GetNtrees());
    for( auto& a : axes ) printf("  %-32s %zu nodes  %g .. %g\n",
                                 a.name.c_str(), a.values.size(), a.values.front(), a.values.back());

    // ---- the scan ---------------------------------------------------------
    TChain gt("FitterEngine/fit/gridTree");
    gt.Add(filePattern.c_str());

    // GridValues/Node replaced Grided/Point; read whichever this file has
    const char* vBranch = gt.GetBranch("GridValues") ? "GridValues" : "Grided";
    double llh, gv[64];
    gt.SetBranchAddress("LLH", &llh);
    gt.SetBranchAddress(vBranch, gv);

    const double BIG = 1e30;
    std::vector<TH1D*> h1(nAx, nullptr);
    std::map<std::pair<int,int>, TH2D*> h2;

    for( int i=0;i<nAx;i++ ){
        h1[i] = new TH1D(Form("prof_%s", axes[i].name.c_str()),
                         Form(";%s;#Delta#chi^{2}", axes[i].label.c_str()),
                         int(axes[i].values.size()), axes[i].edges.data());
        for( int b=1;b<=h1[i]->GetNbinsX();b++ ) h1[i]->SetBinContent(b, BIG);
    }
    for( int i=0;i<nAx;i++ ) for( int j=i+1;j<nAx;j++ ){
        auto* h = new TH2D(Form("prof_%s_%s", axes[i].name.c_str(), axes[j].name.c_str()),
                           Form(";%s;%s", axes[i].label.c_str(), axes[j].label.c_str()),
                           int(axes[i].values.size()), axes[i].edges.data(),
                           int(axes[j].values.size()), axes[j].edges.data());
        for( int bx=1;bx<=h->GetNbinsX();bx++ ) for( int by=1;by<=h->GetNbinsY();by++ )
            h->SetBinContent(bx,by,BIG);
        h2[{i,j}] = h;
    }

    double best = BIG;
    std::vector<double> bestAt(nAx, 0.);
    Long64_t n = gt.GetEntries();
    for( Long64_t e=0;e<n;e++ ){
        gt.GetEntry(e);
        if( llh < best ){ best = llh;
            for( int i=0;i<nAx;i++ ) bestAt[i] = gv[i]*axes[i].scale; }
        for( int i=0;i<nAx;i++ ){
            int b = h1[i]->FindBin( gv[i]*axes[i].scale );
            if( llh < h1[i]->GetBinContent(b) ) h1[i]->SetBinContent(b, llh);
        }
        for( auto& kv : h2 ){
            int i=kv.first.first, j=kv.first.second;
            int b = kv.second->FindBin( gv[i]*axes[i].scale, gv[j]*axes[j].scale );
            if( llh < kv.second->GetBinContent(b) ) kv.second->SetBinContent(b, llh);
        }
    }
    printf("%lld entries, minimum LLH %.5f at", n, best);
    for( int i=0;i<nAx;i++ ) printf("  %s = %g", axes[i].name.c_str(), bestAt[i]);
    printf("\n");

    int empty1=0, empty2=0;
    for( auto* h : h1 ) for( int b=1;b<=h->GetNbinsX();b++ ){
        if( h->GetBinContent(b) >= BIG ){ h->SetBinContent(b,0); empty1++; }
        else h->SetBinContent(b, h->GetBinContent(b)-best); }
    for( auto& kv : h2 ) for( int bx=1;bx<=kv.second->GetNbinsX();bx++ )
        for( int by=1;by<=kv.second->GetNbinsY();by++ ){
            int b = kv.second->GetBin(bx,by);
            if( kv.second->GetBinContent(b) >= BIG ){ kv.second->SetBinContent(b,0); empty2++; }
            else kv.second->SetBinContent(b, kv.second->GetBinContent(b)-best); }
    if( empty1 || empty2 ) printf("unfilled bins: %d in 1D, %d in 2D\n", empty1, empty2);

    // ---- draw -------------------------------------------------------------
    auto out = new TFile(outPath.c_str(), "RECREATE");

    for( int i=0;i<nAx;i++ ){
        auto cv = new TCanvas(Form("c_prof_%s", axes[i].name.c_str()), axes[i].label.c_str(), 900, 650);
        cv->SetLeftMargin(0.13); cv->SetGridx(); cv->SetGridy();
        h1[i]->SetLineColor(kBlack); h1[i]->SetLineWidth(3);
        h1[i]->GetYaxis()->SetTitleSize(0.050); h1[i]->GetYaxis()->SetTitleOffset(1.15);
        h1[i]->GetXaxis()->SetTitleSize(0.050);
        h1[i]->GetYaxis()->SetLabelSize(0.042); h1[i]->GetXaxis()->SetLabelSize(0.042);
        h1[i]->SetMinimum(0);
        h1[i]->Draw("HIST L");

        auto m = new TMarker(bestAt[i], 0., 29);
        m->SetMarkerColor(kRed+1); m->SetMarkerSize(2.0); m->Draw();

        out->cd(); cv->Write(); h1[i]->Write();
    }

    for( auto& kv : h2 ){
        int i=kv.first.first, j=kv.first.second;
        auto* h = kv.second;
        auto cv = new TCanvas(Form("c_%s", h->GetName()), h->GetName(), 900, 700);
        cv->SetLeftMargin(0.13); cv->SetRightMargin(0.16);

        // COLZ skips a bin holding exactly the minimum, which leaves a white
        // hole at the best fit
        h->SetMinimum(-1e-6);
        h->SetMaximum( kContours[kNContours-1] * 2. );
        h->GetZaxis()->SetTitle("#Delta#chi^{2}");
        h->GetYaxis()->SetTitleSize(0.048); h->GetYaxis()->SetTitleOffset(1.35);
        h->GetXaxis()->SetTitleSize(0.048);
        h->GetYaxis()->SetLabelSize(0.040); h->GetXaxis()->SetLabelSize(0.040);
        h->Draw("COLZ");

        // contour levels belong to the drawn object, so the overlay is a clone
        auto hc = (TH2D*) h->Clone(Form("%s_cont", h->GetName()));
        hc->SetContour(kNContours, kContours);
        hc->SetLineColor(kBlack); hc->SetLineWidth(2);
        hc->Draw("CONT3 SAME");

        auto m = new TMarker(bestAt[i], bestAt[j], 29);
        m->SetMarkerColor(kRed+1); m->SetMarkerSize(2.4); m->Draw();

        TString lv;
        for( int k=0;k<kNContours;k++ ) lv += Form("%s%.2f", k?", ":"", kContours[k]);
        auto tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.032);
        tx->DrawLatex(0.15, 0.93, Form("#Delta#chi^{2} contours: %s", lv.Data()));

        out->cd(); cv->Write(); h->Write();
    }

    out->Close();
    printf("wrote %s : %d 1D and %zu 2D profiles\n", outPath.c_str(), nAx, h2.size());
}
