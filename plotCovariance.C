// The prior covariance the fit is given: the four parameter-set blocks placed
// block-diagonally in the order GUNDAM loads them, plus the correlation matrix
// and the prior sigmas. Draws three canvases and writes them with the matrices.
//
//   root -l -b -q "plotCovariance.C(\"$data_dir/parameters\",\"outputs/priorCovariance.root\")"

void plotCovariance(std::string inputs_dir, std::string out_file){

    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetNumberContours(99);
    TH1::AddDirectory(kFALSE);

    struct Block { const char* file; const char* prefix; const char* label; };
    const std::vector<Block> blocks = {
        {"oscProb/oscCovHK.root", "osc",  "Oscillation"},
        {"flux/fluxCovHK.root",   "flux", "Flux"},
        {"det/detCovHK.root",     "det",  "Detector"},
        {"xsec/xsecCovHK.root",   "xsec", "Cross-section"}
    };

    std::vector<TMatrixD> mats;
    std::vector<std::string> names;
    std::vector<int> start;
    int n = 0;
    for( auto& b : blocks ){
        TFile f( (inputs_dir + "/" + b.file).c_str() );
        auto* m  = (TMatrixD*)  f.Get(Form("%s_param_cov",   b.prefix));
        auto* na = (TObjArray*) f.Get(Form("%s_param_names", b.prefix));
        if( !m || !na ){ printf("missing %s_param_cov/names in %s\n", b.prefix, b.file); return; }
        mats.emplace_back(*m);
        start.push_back(n);
        for( int i=0; i<na->GetEntries(); i++ ) names.push_back(na->At(i)->GetName());
        n += m->GetNrows();
    }
    printf("%d parameters:", n);
    for( size_t k=0; k<blocks.size(); k++ ) printf("  %s %d", blocks[k].label, mats[k].GetNrows());
    printf("\n");

    // block-diagonal total, and its correlation
    TMatrixD cov(n,n), cor(n,n);
    for( size_t k=0; k<mats.size(); k++ )
        for( int i=0; i<mats[k].GetNrows(); i++ )
            for( int j=0; j<mats[k].GetNcols(); j++ )
                cov(start[k]+i, start[k]+j) = mats[k](i,j);
    for( int i=0; i<n; i++ )
        for( int j=0; j<n; j++ )
            cor(i,j) = cov(i,j) / std::sqrt(cov(i,i)*cov(j,j));

    auto hCov = new TH2D("priorCovariance", ";;", n,0,n, n,0,n);
    auto hCor = new TH2D("priorCorrelation",";;", n,0,n, n,0,n);
    auto hSig = new TH1D("priorSigma", ";;prior #sigma", n,0,n);
    for( int i=0; i<n; i++ ){
        hCov->GetXaxis()->SetBinLabel(i+1, names[i].c_str());
        hCov->GetYaxis()->SetBinLabel(i+1, names[i].c_str());
        hCor->GetXaxis()->SetBinLabel(i+1, names[i].c_str());
        hCor->GetYaxis()->SetBinLabel(i+1, names[i].c_str());
        hSig->GetXaxis()->SetBinLabel(i+1, names[i].c_str());
        hSig->SetBinContent(i+1, std::sqrt(cov(i,i)));
        for( int j=0; j<n; j++ ){
            hCov->SetBinContent(i+1, j+1, cov(i,j));
            hCor->SetBinContent(i+1, j+1, cor(i,j));
        }
    }

    auto dress = [&](TH2D* h){
        h->GetXaxis()->SetLabelSize(0.012); h->GetYaxis()->SetLabelSize(0.012);
        h->GetXaxis()->LabelsOption("v");
        h->GetZaxis()->SetLabelSize(0.02);
        // block boundaries and labels
        for( size_t k=1; k<start.size(); k++ ){
            auto lx = new TLine(start[k],0,start[k],n); lx->SetLineColor(kBlack); lx->SetLineWidth(2); lx->Draw();
            auto ly = new TLine(0,start[k],n,start[k]); ly->SetLineColor(kBlack); ly->SetLineWidth(2); ly->Draw();
        }
        // narrow blocks get their label on a second row so neighbours do not collide
        for( size_t k=0; k<blocks.size(); k++ ){
            double c = start[k] + 0.5*mats[k].GetNrows();
            double y = n + ( mats[k].GetNrows() < 10 ? 2.6 : 0.8 );
            auto t = new TLatex(c, y, blocks[k].label);
            t->SetTextAlign(21); t->SetTextSize(0.022); t->Draw();
        }
    };

    auto out = new TFile(out_file.c_str(), "RECREATE");

    auto c1 = new TCanvas("c_priorCorrelation","prior correlation",1400,1300);
    c1->SetLeftMargin(0.17); c1->SetBottomMargin(0.17); c1->SetRightMargin(0.12); c1->SetTopMargin(0.05);
    auto exCor = new TExec("exCor","gStyle->SetPalette(kTemperatureMap);");
    hCor->SetMinimum(-1); hCor->SetMaximum(1);
    hCor->Draw("AXIS"); exCor->Draw(); hCor->Draw("COLZ SAME"); dress(hCor);
    out->cd(); c1->Write(); hCor->Write();

    // entries span 1e-9 (dm2) to 1.4 (dCP), so log z; ROOT leaves non-positive
    // entries blank on a log scale, and the correlation canvas carries the sign
    auto c2 = new TCanvas("c_priorCovariance","prior covariance",1400,1300);
    c2->SetLeftMargin(0.17); c2->SetBottomMargin(0.17); c2->SetRightMargin(0.12); c2->SetTopMargin(0.05);
    c2->SetLogz();
    double zmin = 1e300;
    for( int i=0; i<n; i++ ) if( cov(i,i) > 0 && cov(i,i) < zmin ) zmin = cov(i,i);
    hCov->SetMinimum(zmin);
    auto exCov = new TExec("exCov","gStyle->SetPalette(kBird);");
    hCov->Draw("AXIS"); exCov->Draw(); hCov->Draw("COLZ SAME"); dress(hCov);
    out->cd(); c2->Write(); hCov->Write();

    auto c3 = new TCanvas("c_priorSigma","prior sigma",1400,600);
    c3->SetLogy(); c3->SetBottomMargin(0.30); c3->SetLeftMargin(0.08); c3->SetGridy();
    hSig->GetXaxis()->SetLabelSize(0.018); hSig->GetXaxis()->LabelsOption("v");
    hSig->GetYaxis()->SetTitleSize(0.045); hSig->GetYaxis()->SetTitleOffset(0.8);
    hSig->SetFillColor(kAzure-9); hSig->SetLineColor(kAzure+2);
    hSig->Draw("HIST");
    for( size_t k=1; k<start.size(); k++ ){
        auto l = new TLine(start[k], hSig->GetMinimum(), start[k], hSig->GetMaximum()*1.5);
        l->SetLineColor(kBlack); l->SetLineWidth(2); l->Draw();
    }
    out->cd(); c3->Write(); hSig->Write();

    cov.Write("priorCovarianceMatrix");
    cor.Write("priorCorrelationMatrix");
    out->Close();
    printf("wrote %s\n", out_file.c_str());
}
