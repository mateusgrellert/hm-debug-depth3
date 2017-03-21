#include "TComDebug.h"


std::string TComDebug::debugFilePath;
FILE* TComDebug::debugFile;
std::map<std::string, std::vector<double> > TComDebug::statsMap; 
std::vector<std::string > TComDebug::cuOrderMap; 
bool TComDebug::encodeStarted;
int TComDebug::numCodedPics;
double TComDebug::ctx_ctu_depth;
double TComDebug::lambda;
double TComDebug::SAD;
double TComDebug::MSE;
double TComDebug::SSE;
double TComDebug::nonZeroCoeff;
double TComDebug::cost_2Nx2N;
double TComDebug::cost_2NxN;
double TComDebug::cost_Nx2N;
double TComDebug::cost_MSM;
double TComDebug::delta_QP;
double TComDebug::cu_time;

void TComDebug::getMSE( TComYuv *pcYuvSrc0, TComYuv *pcYuvSrc1,int width){
    const ComponentID compID=ComponentID(COMPONENT_Y);
    UInt uiTrUnitIdx = 0;
    MSE = 0;
    SAD = 0;
    SSE = 0;
    int diff, power2diff;
    const Int uiPartWidth = width >> (pcYuvSrc0->getComponentScaleX(compID));
    const Int uiPartHeight= width >> (pcYuvSrc0->getComponentScaleY(compID));

    Pel* pSrc0 = pcYuvSrc0->getAddr( compID, uiTrUnitIdx, uiPartWidth );
    Pel* pSrc1 = pcYuvSrc1->getAddr( compID, uiTrUnitIdx, uiPartWidth );
          

    Int  iSrc0Stride = pcYuvSrc0->getStride(compID);
    Int  iSrc1Stride = pcYuvSrc1->getStride(compID);

    for (Int y = uiPartHeight-1; y >= 0; y-- )
    {
      for (Int x = uiPartWidth-1; x >= 0; x-- )
      {
        diff = abs(pSrc0[x] - pSrc1[x]);
        power2diff = diff*diff;
        SAD += diff;
        SSE += power2diff;
        MSE += power2diff;
      }
      pSrc0 += iSrc0Stride;
      pSrc1 += iSrc1Stride;
    }
    MSE = MSE/(uiPartHeight*uiPartWidth);
}

int TComDebug::getFME(TComDataCU *&cu){
    
        int fme = 0;
    int mv_h = 0;
    int mv_v = 0;
    double frac_mv_h = 0;
    double frac_mv_v = 0;
    


        int uiRefListIdx = 0;

    for ( uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
        {
            if ( cu->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
                mv_h = cu->getCUMvField(RefPicList(uiRefListIdx))->getMv( 0 ).getHor();
                mv_v = cu->getCUMvField(RefPicList(uiRefListIdx))->getMv( 0 ).getVer();
                frac_mv_h = (abs(mv_h) % 4)/4.0;
                frac_mv_v = (abs(mv_v) % 4)/4.0;

                if ((frac_mv_v == 0.25 ) || (frac_mv_h == 0.25)){
                    fme = 2;
                }                
                else if ((frac_mv_v == 0.75) || (frac_mv_h == 0.75)){
                    fme = 2;
                }
                else if ((frac_mv_v == 0.5) || (frac_mv_h == 0.5)){
                    fme = 1;
                }
                else
                    fme = 0;
                break;
            }
        }
        return fme;
}

double TComDebug::getAverageNeighborDepth(TComDataCU *&cu){
    int numValidCtu = 0;
    TComDataCU* upCtu = cu->getCtuAbove();
    TComDataCU* leftCtu = cu->getCtuLeft();
    TComDataCU* upLeftCtu = cu->getCtuAboveLeft();
    TComDataCU* upRightCtu = cu->getCtuAboveRight();
    TComDataCU* colocCtuL0 = cu->getCUColocated(REF_PIC_LIST_0);
    TComDataCU* colocCtuL1 = cu->getCUColocated(REF_PIC_LIST_1);
    if (colocCtuL0 == colocCtuL1)
        colocCtuL1 = NULL;

    double avgCtxDepth = 0;

    if(upCtu){
        numValidCtu++;
        avgCtxDepth += calcAverageDepthIterative(upCtu);
    }
    if(leftCtu){
        numValidCtu++;
        avgCtxDepth += calcAverageDepthIterative(leftCtu);
    }
    if(upLeftCtu){
        numValidCtu++;
        avgCtxDepth += calcAverageDepthIterative(upLeftCtu);
    }
    if(upRightCtu){
        numValidCtu++;
        avgCtxDepth += calcAverageDepthIterative(upRightCtu);
    }
    if(colocCtuL0){
        numValidCtu++;
        avgCtxDepth += calcAverageDepthIterative(colocCtuL0);
    }
    if(colocCtuL1){
        numValidCtu++;
        avgCtxDepth += calcAverageDepthIterative(colocCtuL1);
    }
    return avgCtxDepth/numValidCtu;
}    
    
double TComDebug::calcAverageDepthIterative(TComDataCU *&cu){
    int uiDepth,uiAbsPartIdx0,uiAbsPartIdx1,uiAbsPartIdx2;
    UInt uiQNumParts0 ;
    UInt uiQNumParts1 ;
    UInt uiQNumParts2 ;
    UInt uiPartUnitIdx0,uiPartUnitIdx1,uiPartUnitIdx2;
    double avgDepth = 0;

    uiAbsPartIdx0 = 0;
    uiAbsPartIdx1 = 0;
    uiAbsPartIdx2 = 0;
    uiDepth = 0;
    uiQNumParts0 = ( cu->getPic()->getNumPartitionsInCtu() >> (uiDepth<<1) )>>2;
    int totalCus = 0;
    for ( uiPartUnitIdx0 = 0; uiPartUnitIdx0 < 4; uiPartUnitIdx0++, uiAbsPartIdx0+=uiQNumParts0 ){
        if (cu->getDepth(uiAbsPartIdx0) == uiDepth){
            avgDepth += uiDepth;
            totalCus ++;
        }
        else{
            uiDepth = 1;
            uiQNumParts1 = ( cu->getPic()->getNumPartitionsInCtu() >> (uiDepth<<1) )>>2;
            for ( uiAbsPartIdx1 = uiAbsPartIdx0, uiPartUnitIdx1 = 0; uiPartUnitIdx1 < 4; uiPartUnitIdx1++, uiAbsPartIdx1+=uiQNumParts1 ){
                if (cu->getDepth(uiAbsPartIdx1) == uiDepth){
                    avgDepth += uiDepth;
                    totalCus ++;
                }
                else{
                    uiDepth = 2;
                    uiQNumParts2 = ( cu->getPic()->getNumPartitionsInCtu() >> (uiDepth<<1) )>>2;
                    for ( uiAbsPartIdx2 = uiAbsPartIdx1, uiPartUnitIdx2 = 0; uiPartUnitIdx2 < 4; uiPartUnitIdx2++, uiAbsPartIdx2+=uiQNumParts2 ){
                        if (cu->getDepth(uiAbsPartIdx2) == uiDepth){
                            avgDepth += uiDepth;
                            totalCus ++;
                        }
                        else{
                            uiDepth = 3;
                            if (cu->getDepth(uiAbsPartIdx2) == uiDepth){
                                avgDepth += uiDepth;
                                totalCus ++;
                            }
                            else{
                                printf("Error!\n");
                                exit(1);
                            }
                        }
                    }
                }
            }
        }
    }
    
    return avgDepth/totalCus;

}
double TComDebug::getFeatureValue(TComDataCU *&cu, int feat_idx){
    double feat_val = 0;
    
    int cmode;
    int cux = cu->getCUPelX();
    int cuy = cu->getCUPelY();
    
    if(cu->isSkipped(0)) cmode = 0;
    else if (cu->isInter(0) && !cu->getQtRootCbf( 0 )) cmode = 1;
    else if (cu->isInter(0)) cmode = 2;
    else cmode = 3;

    
    int colocSplit = 0;
    double colocDepth = 0.0;
    int depth = cu->getDepth(0);
    int width = 64 >> depth;
    
    double ratio2Nx2N_MSM = 0,relRatio2Nx2N_MSM = 0,ratioBest_2Nx2N = 0,ratioBest_MSM = 0,relRatioBest_MSM = 0;
    
    double costBest = cu->getTotalCost();
    if (cost_MSM == 0){
        cost_MSM = 0.1;
    }
    if (cost_2Nx2N == 0)
        cost_2Nx2N = 0.1;
    if (cost_2NxN == 0)
        cost_2NxN = 0.1;
    if (cost_Nx2N == 0)
        cost_Nx2N = 0.1;
    
    if(feat_idx >= 20 && feat_idx <= 24){
        ratio2Nx2N_MSM = cost_2Nx2N / cost_MSM;
        relRatio2Nx2N_MSM = abs((cost_2Nx2N - cost_MSM) / cost_MSM);
        ratioBest_MSM = costBest / cost_MSM;
        relRatioBest_MSM =  abs((costBest - cost_MSM) / cost_MSM); 
        ratioBest_2Nx2N = costBest / cost_2Nx2N;
    }

    int fme = 0;
    int mv_h = 0;
    int mv_v = 0;
    double frac_mv_h = 0;
    double frac_mv_v = 0;
    
    int mvd_h = 0;
    int mvd_v = 0;
    double frac_mvd_h = 0;
    double frac_mvd_v = 0;
    
    int pred_mv_h = 0;
    int pred_mv_v = 0;
    double frac_pred_mv_h = 0;
    double frac_pred_mv_v = 0;
    
    double int_mv_mod = 0;
    double int_pred_mv_mod = 0;
    double int_mvd_mod = 0;
    
    double frac_mv_mod = 0;
    double frac_mvd_mod = 0;
    double frac_pred_mv_mod = 0;
    
    
        
    TComDataCU *colocCuL0 = cu->getCUColocated(REF_PIC_LIST_0);
    TComDataCU *colocCuL1 = cu->getCUColocated(REF_PIC_LIST_1);
    
    if(colocCuL0 == colocCuL1)
        colocCuL1 = NULL;
    
    int n_coloc = 0;
    UInt partIdx = cu->getZorderIdxInCtu();
    if (colocCuL0){
        colocSplit += (int) ((int) (colocCuL0->getDepth(partIdx)) > (int)  (cu->getDepth(0)));
        colocDepth += (int) (colocCuL0->getDepth(partIdx));
        n_coloc++;
    }
    if (colocCuL1){
        colocSplit += (int) ((int) (colocCuL1->getDepth(partIdx)) > (int)  (cu->getDepth(0)));
        colocDepth += (int) (colocCuL1->getDepth(partIdx));
        n_coloc++;
    }
    colocDepth = (double) (colocDepth)/n_coloc;
    
    
    int upSplt =  cu->getAboveSplitFlag( 0, depth );
    int leftSplt =  cu->getLeftSplitFlag( 0, depth ); 
    int upLeftSplt =  cu->getUpLeftSplitFlag( 0, depth ); 
    int upRightSplt =  cu->getUpRightSplitFlag( 0, depth ); 
    int ctxSplit = upSplt + leftSplt + upLeftSplt + upRightSplt + colocSplit;
    
    int upDepth =  cu->getAboveDepth( 0, depth );
    int leftDepth =  cu->getLeftDepth( 0, depth ); 
    int upLeftDepth =  cu->getUpLeftDepth( 0, depth ); 
    int upRightDepth =  cu->getUpRightDepth( 0, depth ); 
    double ctxDepth = (upDepth + leftDepth + upLeftDepth + upRightDepth + colocDepth)/(4+n_coloc); 
    
    int uiRefListIdx = 0;
    
    if(feat_idx >= 31 && feat_idx <= 39){
        for ( uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
        {
            if ( cu->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
                mv_h = cu->getCUMvField(RefPicList(uiRefListIdx))->getMv( 0 ).getHor();
                mv_v = cu->getCUMvField(RefPicList(uiRefListIdx))->getMv( 0 ).getVer();
                frac_mv_h = (abs(mv_h) % 4)/4.0;
                frac_mv_v = (abs(mv_v) % 4)/4.0;

                mvd_h = cu->getCUMvField(RefPicList(uiRefListIdx))->getMvd( 0 ).getHor();
                mvd_v = cu->getCUMvField(RefPicList(uiRefListIdx))->getMvd( 0 ).getVer();
                frac_mvd_h = (abs(mvd_h) % 4)/4.0;
                frac_mvd_v = (abs(mvd_v) % 4)/4.0;

                pred_mv_h = mvd_h + mv_h;
                pred_mv_v = mvd_v + mv_v;
                frac_pred_mv_h = (abs(pred_mv_h) % 4)/4.0;
                frac_pred_mv_v = (abs(pred_mv_v) % 4)/4.0;

                int_mv_mod = sqrt(((mv_h >> 2) * (mv_h >> 2)) + ((mv_v >> 2) * (mv_v >> 2)));
                int_pred_mv_mod = sqrt(((pred_mv_h >> 2) * (pred_mv_h >> 2)) + ((pred_mv_v >> 2) * (pred_mv_v >> 2)));
                int_mvd_mod = sqrt(((mvd_h >> 2) * (mvd_h >> 2)) + ((mvd_v >> 2) * (mvd_v >> 2)));

                frac_mv_mod = sqrt((frac_mv_h * frac_mv_h) + (frac_mv_v * frac_mv_v));
                frac_mvd_mod = sqrt((frac_mvd_h * frac_mvd_h) + (frac_mvd_v * frac_mvd_v));
                frac_pred_mv_mod = sqrt((frac_pred_mv_h * frac_pred_mv_h) + (frac_pred_mv_v * frac_pred_mv_v));

                if ((frac_mv_v == 0.25) || (frac_mv_h == 0.25)){
                    fme = 2;
                }                
                else if ((frac_mv_v == 0.75) || (frac_mv_h == 0.75)){
                    fme = 2;
                }
                else if ((frac_mv_v == 0.5) || (frac_mv_h == 0.5)){
                    fme = 1;
                }
                else
                    fme = 0;
                break;
            }
        }
    }
    

    
    switch(feat_idx){
        case 0: feat_val = cu->getPic()->getPOC(); break;
        case 1: feat_val = cux; break;
        case 2: feat_val = cuy; break;
        case 3: feat_val =  width; break;
        case 4: feat_val =  cu->getDepth(0); break;
        case 5: feat_val =  delta_QP; break;
        case 6: feat_val = cmode; break;
        case 7: feat_val = cu->getPartitionSize(0); break;
        case 8: feat_val = cu->getTotalBits(); break;
        case 9: feat_val = cu->getTotalDistortion(); break;
        case 10: feat_val = cu->getTotalCost(); break;
        case 11: feat_val = cu->getTotalBins(); break;
        case 12: feat_val  = cu->getTotalBits()/ cu->getTotalCost() * lambda; break ;// lambda
        
        case 13: feat_val = SAD; break;
        case 14: feat_val = SSE; break;
        case 15: feat_val = MSE; break;
        case 16: feat_val =  cost_2Nx2N; break;
        case 17: feat_val =  cost_MSM; break;
        case 18: feat_val =  cost_2NxN; break;
        case 19: feat_val =  cost_Nx2N; break;
        case 20: feat_val =  ratio2Nx2N_MSM; break;
        case 21: feat_val =  relRatio2Nx2N_MSM; break;
        case 22: feat_val =  ratioBest_2Nx2N; break;
        case 23: feat_val =  ratioBest_MSM; break;
        case 24: feat_val =  relRatioBest_MSM; break;
        case 25: feat_val =  cu->getTransformIdx(0); break;
        case 26: feat_val =  cu->getQtRootCbf( 0 );  break;
        case 27: feat_val =  cu->getTransformSkip(0,COMPONENT_Y); break;
        case 28: feat_val =  nonZeroCoeff; break;
        case 29: feat_val = cu_time * 1000.0; break;
        case 30: feat_val = cu->getCUMvField(RefPicList(uiRefListIdx))->getRefIdx(0) ; break;
        case 31: feat_val = int_mv_mod; break;
        case 32: feat_val = frac_mv_mod; break;
        case 33: feat_val = int_pred_mv_mod; break;
        case 34: feat_val = frac_pred_mv_mod; break;
        case 35: feat_val = int_mvd_mod; break;
        case 36: feat_val = frac_mvd_mod; break;
        case 37: feat_val = cu->getMVPIdx(RefPicList( uiRefListIdx) , 0); break;
        case 38: feat_val = cu->getInterDir(0); break;
        case 39: feat_val = fme; break;
        case 40: feat_val = colocSplit; break;
        case 41: feat_val = upSplt; break;
        case 42: feat_val = leftSplt; break;
        case 43: feat_val = upLeftSplt; break;
        case 44: feat_val = upRightSplt; break;
        case 45: feat_val = ctxSplit; break;
                
        case 46: feat_val = colocDepth; break;
        case 47: feat_val = upDepth; break;
        case 48: feat_val = leftDepth; break;
        case 49: feat_val = upLeftDepth; break;
        case 50: feat_val = upRightDepth; break;
        case 51: feat_val = ctxDepth; break;
        
        case 52: feat_val = ctx_ctu_depth; break;
           
        default: fprintf(stderr,"Error! Feature IDX %d not supported\n", feat_idx); exit(1) ;
        }
    return feat_val;
}

void TComDebug::EVAL_CU(TComDataCU *&cu){
    std::string sstr_str = TComDebug::getMapString(cu, cu->getDepth(0), 0);

    cuOrderMap.push_back(sstr_str);
    int feat_idx;
    for(feat_idx = 0; feat_idx < N_FEATURES-N_LABELS; feat_idx++){
        statsMap[sstr_str].push_back(getFeatureValue(cu, feat_idx));
    }
    
  
}


std::string TComDebug::getMapString(TComDataCU *cu, int depth, int partIdx){
    std::stringstream sstr;
    
    int cu_x = cu->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[partIdx] ];
    int cu_y = cu->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[partIdx] ];
    int width = 64 >> depth;

    //int cu_x = pcCU->getCUPelX() ;
    //int cu_y = pcCU->getCUPelY() ;
    sstr << cu->getPic()->getPOC() << "_" << cu_x << "x" << cu_y << "_" << width;
    
    return sstr.str();
}


void TComDebug::writeStats(){
    //std::map<std::string, std::vector<double> >::iterator it_i;
    //std::vector<double>::iterator it_j;
    
    int i,j;
    for(i = 0; i < cuOrderMap.size(); i++){
        const std::string &s = cuOrderMap[i];
        if (statsMap[s].size() < N_FEATURES) // skipping CUs below the ones that were not split
            continue;
        for(j = 0; j < statsMap[s].size()-1; j++ )
            fprintf(debugFile,"%f,",(double) statsMap[s][j]);
        fprintf(debugFile,"%f\n",(double) statsMap[s][j]);
    }
}

void TComDebug::printHeader(){

        int i;

    for(i = 0; i < sizeof(featureMap)/sizeof(featureMap[0])-1; i++)
        fprintf(debugFile,"%s,",featureMap[i].c_str());
    fprintf(debugFile,"%s\n",featureMap[i].c_str());
}