#ifndef QCDNUM_H
#define QCDNUM_H

#include <string>
#include <iostream>
using namespace std;

namespace MBUTIL {

    /************************************************************/
    /*  MBUTIL comparison routines from compa.f                 */
    /************************************************************/

    int lmb_eq(double a, double b, double epsi);
    int lmb_ne(double a, double b, double epsi);
    int lmb_ge(double a, double b, double epsi);
    int lmb_le(double a, double b, double epsi);
    int lmb_gt(double a, double b, double epsi);
    int lmb_lt(double a, double b, double epsi);

    /************************************************************/
    /*  MBUTIL hash routines from hash.f                        */
    /************************************************************/

    int imb_ihash(int iseed, int *imsg, int n);
    int imb_jhash(int iseed, double *dmsg, int n);
    int imb_dhash(int iseed, double *dmsg, int n);

    /************************************************************/
    /*  MBUTIL character routines from mchar.f                  */
    /************************************************************/

    void smb_itoch(int in, string &chout, int &leng);
    void smb_dtoch(double dd, int n, string &chout, int &leng);
    void smb_hcode(int ihash, string &hcode);

    /************************************************************/
    /*  MBUTIL interpolation routines from polint.f             */
    /************************************************************/

    // Already defined in file wspace.f
    // inline int iaFtoC(int ia) { return ia-1; };
    // inline int iaCtoF(int ia) { return ia+1; };

    void smb_spline(int n, double *x, double *y, double *b, double *c, double *d);
    double dmb_seval(int n, double u, double *x, double *y, double *b, double *c, double *d);

    /************************************************************/
    /*  MBUTIL utility routines from utils.f                    */
    /************************************************************/

    int imb_version();
    double dmb_gauss(double (*f)(double*), double a, double b,
                     double &eps);
    double dmb_gaus1(double (*f)(double*), double a, double b);
    double dmb_gaus2(double (*f)(double*), double a, double b);
    double dmb_gaus3(double (*f)(double*), double a, double b);
    double dmb_gaus4(double (*f)(double*), double a, double b);

    /************************************************************/
    /*  MBUTIL vector routines from vectors.f                   */
    /************************************************************/

    void smb_ifill(int *ia, int n, int ival);
    void smb_vfill(double *a, int n, double val);
    void smb_vmult(double *a, int n, double val);
    void smb_vcopy(double *a, double* b, int n);
    void smb_icopy(int *ia, int* ib, int n);
    void smb_vitod(int *ia, double *b, int n);
    void smb_vdtoi(double *a, int *ib, int n);
    void smb_vaddv(double *a, double *b, double *c, int n);
    void smb_vminv(double *a, double *b, double *c, int n);
    double dmb_vdotv(double *a, double *b, int n);
    double dmb_vnorm(int m, double *a, int n);
    int lmb_vcomp(double *a, double *b, int n, double epsi);
    double dmb_vpsum(double *a, double *w, int n);
}

namespace WSTORE {

    /************************************************************/
    /*  WSTORE istore routines from istore.f                    */
    /************************************************************/

    // alrady defined in file wstore.f
    // inline int iaFtoC(int ia) { return ia-1; };
    // inline int iaCtoF(int ia) { return ia+1; };

    void sws_IwInit(int *iw, int nw, string txt);
    void sws_SetIwN(int *iw, int nw);
    int iws_Iarray(int* iw, int imi, int ima);
    int iws_IAread(int* iw, int* iarr, int n);
    int iws_DAread(int* iw, double* darr, int n);
    void sws_IwWipe(int* obj);
    int iws_IhSize();

    int iws_IsaIstore(int *obj);
    int iws_IwSize(int *obj);
    int iws_IwNused(int *obj);
    int iws_IwNarrays(int *obj);
    int iws_IaLastObj(int *obj);
    int iws_IwNheader(int *obj);
    int iws_IwObjectType(int *obj);
    int iws_IwAsize(int *obj);
    int iws_IwAnumber(int *obj);
    int iws_IwAfprint(int *obj);
    int iws_IwAdim(int *obj);
    int iws_IwAimin(int *obj);
    int iws_IwAimax(int *obj);
    int iws_IaAbegin(int *obj);
    int iws_IaAend(int *obj);
    int iws_IwKnul(int *obj);
    int iws_ArrayI(int *obj, int i);
    void sws_IwTree(int *iw);
    void sws_IwHead(int *obj);

    /************************************************************/
    /*  WSTORE workspace routines from wstore.f                 */
    /************************************************************/

    inline int iaFtoC(int ia) { return ia-1; };
    inline int iaCtoF(int ia) { return ia+1; };

    int iws_Version();
    int iws_WsInit(double *w, int nw, int nt, string txt);
    void sws_SetWsN(double *w, int nw);
    void sws_Stampit(double *w);
    int iws_NewSet(double *w);
    int iws_WTable(double* w, int* imi, int* ima, int n);
    int iws_WClone(double *w1, int ia1, double *w2);
    void sws_TbCopy(double *w1, int ia1,
                    double *w2, int ia2, int itag);
    void sws_WsWipe(double *w, int ia);
    void sws_TsDump(string fname, int key,
                    double *w, int ia, int &ierr);
    int iws_TsRead(string fname, int key, double *w, int &ierr);
    void sws_WsMark(int &mws, int &mset, int &mtab );
    int iws_HdSize();
    int iws_TbSize(int *imi, int *ima, int ndim);
    int iws_Tpoint(double *w, int ia, int* index, int n);

    int iws_IsaWorkspace(double *w);
    int iws_SizeOfW(double *w);
    int iws_WordsUsed(double *w);
    int iws_Nheader(double *w);
    int iws_Ntags(double *w);
    int iws_HeadSkip(double *w);
    int iws_IaRoot();
    int iws_IaDrain(double *w);
    int iws_IaNull(double *w);
    int iws_ObjectType(double *w, int ia);
    int iws_ObjectSize(double *w, int ia);
    int iws_Nobjects(double *w, int ia);
    int iws_ObjectNumber(double *w, int ia);
    int iws_FingerPrint(double *w, int ia);
    int iws_IaFirstTag(double *w, int ia);
    int iws_TableDim(double *w, int ia);
    int iws_IaKARRAY(double *w, int ia);
    int iws_IaIMIN(double *w, int ia);
    int iws_IaIMAX(double *w, int ia);
    int iws_BeginTbody(double *w, int ia);
    int iws_EndTbody(double *w, int ia);
    int iws_TFskip(double *w, int ia);
    int iws_TBskip(double *w, int ia);
    int iws_SFskip(double *w, int ia);
    int iws_SBskip(double *w, int ia);
    void sws_WsTree(double *w);
    void sws_WsHead(double *w, int ia);
}

namespace SPLINT {


    /************************************************************/
    /*  SPLINT routines from usrsplint.f                        */
    /************************************************************/

    int isp_spvers();
    void ssp_spinit(int nuser);
    void ssp_uwrite(int i, double val);
    double dsp_uread(int i);
    int isp_sxmake(int istep);
    int isp_sxuser(double* xarr, int nx);
    void ssp_sxfill(int ias, double (*fun)(int*,int*,bool*),
                    int iq);
    void ssp_sxfpdf(int ias, int iset, double *def, int isel,
                    int iq);
    void ssp_sxf123(int ias, int iset, double *def, int istf,
                    int iq);
    int isp_sqmake(int istep);
    int isp_squser(double* qarr, int nq);
    void ssp_sqfill(int ias, double (*fun)(int*,int*,bool*),
                    int ix);
    void ssp_sqfpdf(int ias, int iset, double *def, int isel,
                    int ix);
    void ssp_sqf123(int ias, int iset, double *def, int istf,
                    int ix);
    double dsp_funs1(int ia, double z, int ichk);
    double dsp_ints1(int ia, double u1, double u2);
    int isp_s2make(int istepx, int istepq);
    int isp_s2user(double *xarr, int nx, double *qarr, int nq);
    void ssp_s2fill(int ias, double (*fun)(int*,int*,bool*),
                    double rs);
    void ssp_s2fpdf(int ias, int iset, double *def, int isel,
                    double rs);
    void ssp_s2f123(int ias, int iset, double *def, int istf,
                    double rs);
    double dsp_funs2(int ia, double x, double q, int ichk);
    double dsp_ints2(int ia, double x1, double x2, double q1,
                     double q2, double rs, int np);
    void ssp_extrapu(int ia, int n);
    void ssp_extrapv(int ia, int n);
    int isp_splinetype(int ia);
    void ssp_splims(int ia, int &nu, double &umi, double &uma,
                    int &nv, double &vmi, double &vma, int &nb);
    void ssp_unodes(int ia, double* array, int n, int &nus);
    void ssp_vnodes(int ia, double* array, int n, int &nvs);
    void ssp_nprint(int ia);
    double dsp_rscut(int ia);
    double dsp_rsmax(int ia, double rs);
    void ssp_mprint();
    int isp_spsize(int ia);
    void ssp_erase(int ia);
    void ssp_spdump(int ia, string fname);
    int isp_spread(string fname);
    void ssp_spsetval(int ia, int i, double val);
    double dsp_spgetval(int ia, int i);
}

namespace QCDNUM {


    /************************************************************/
    /*  QCDNUM address function from usrGlobalId.f              */
    /************************************************************/

    int ipdftab(int iset, int id);

    /************************************************************/
    /*  QCDNUM convolution engine from usrcvol.f                */
    /************************************************************/

    void makewta(double *w, int jd,
                 double (*afun)(double*,double*,int*),
                 double (*achi)(double*));
    void makewtb(double *w, int jd,
                 double (*bfun)(double*,double*,int*),
                 double (*achi)(double*),int ndel);
    void makewrs(double *w, int jd,
                 double (*rfun)(double*,double*,int*),
                 double (*sfun)(double*,double*,int*),
                 double (*achi)(double*),int ndel);
    void makewtd(double *w, int jd,
                 double (*dfun)(double*,double*,int*),
                 double (*achi)(double*));
    void makewtx(double *w, int jd);
    int idspfun(string pname, int iord, int jset);
    void scalewt(double *w, double c, int jd);
    void copywgt(double *w, int jd1, int jd2, int iadd);
    void wcrossw(double *w, int jda, int jdb, int jdc, int iadd);
    void wtimesf(double *w, double (*fun)(int*,int*),
                 int jd1, int jd2, int iadd);
    double fcrossk(double *w, int idw, int jdum, int idf,
                   int ix, int iq);
    double fcrossf(double *w, int idw, int jdum, int ida, int idb,
                   int ix, int iq);
    void efromqq(double *qvec, double *evec, int nf);
    void qqfrome(double *evec, double *qvec, int nf);
    void stfunxq(double (*fun)(int*,int*), double *x, double *q,
                 double *f, int n, int jchk);

    /************************************************************/
    /*  QCDNUM error message routines from usrerr.f             */
    /************************************************************/

    void setumsg(string name);
    void clrumsg();

    /************************************************************/
    /*  QCDNUM toolbox evolution from usrevnn.f                 */
    /************************************************************/

    void evfilla(double *w, int id, double (*func)(int*,int*,int*));
    double evgetaa(double *w, int id, int iq, int &nf, int &ithresh);
    void evdglap(double *w, int idw, int ida, int idf, double *start, int m, int n, int *iqlim, int &nf, double &eps);
    void cpyparw(double *w, double *array, int n, int kset);
    void useparw(double *w, int kset);
    int keyparw(double *w, int kset);
    int keygrpw(double *w, int kset, int igroup);
    double evpdfij(double *w, int id, int ix, int jq, int jchk);
    void evplist(double *w, int id, double *x, double *q, double *f, int n, int jchk);
    void evtable(double *w, int id, double *xx, int nx, double *qq, int nq, double *pdf, int jchk);
    void evpcopy(double *w, int id, double *def, int n, int jset);

    /************************************************************/
    /*  QCDNUM evolution routines from usrevol.f                */
    /************************************************************/

    double rfromf(double fscale2);
    double ffromr(double rscale2);
    double asfunc(double r2, int &nf, int &ierr);
    double altabn(int iset, int iq, int n, int &ierr);
    void evolfg(int iset, double (*func)(int*, double*), double *def, int iq0, double &epsi);
    void evsgns(int iset, double (*func)(int*, double*), int *isns, int n, int iq0, double &epsi);
    void pdfcpy(int iset1, int iset2);
    void extpdf(double (*func)(int*, double*, double*, bool*), int iset, int n, double offset, double &epsi);
    void usrpdf(double (*func)(int*, double*, double*, bool*), int iset, int n, double offset, double &epsi);
    int nptabs(int iset);
    int ievtyp(int iset);

    /************************************************************/
    /*  QCDNUM fast convolution routines from usrfast.          */
    /************************************************************/

    void idscope(double *w, int jset);
    void fastini(double *xlist, double *qlist, int n, int jchk);
    void fastclr(int ibuf);
    void fastinp(double *w, int jdf, double *coef, int jbuf, int iadd);
    void fastepm(int jdum, int jdf, int jbuf);
    void fastsns(int jset, double *qvec, int isel, int jbuf);
    void fastsum(int jset, double *coef, int jbuf);
    void fastfxk(double *w, int *idwt, int ibuf1, int jbuf2);
    void fastfxf(double *w, int idx, int ibuf1, int ibuf2, int jbuf3);
    void fastkin(int ibuf, double (*fun)(int*,int*,int*,int*));
    void fastcpy(int ibuf1, int ibuf2, int iadd);
    void fastfxq(int ibuf, double *stf, int n);

    /************************************************************/
    /*  QCDNUM grid routines from usrgrd.f                      */
    /************************************************************/

    void gxmake(double *xmin, int *iwt, int n, int nxin, int &nxout, int iord);
    int ixfrmx(double x);
    int xxatix(double x, int ix);
    double xfrmix(int ix);
    void gxcopy(double *array, int n, int &nx);
    void gqmake(double *qarr, double *wgt, int n, int nqin, int &nqout);
    int iqfrmq(double q2);
    int qqatiq(double q2, int iq);
    double qfrmiq(int iq);
    void gqcopy(double *array, int n, int &nq);
    void grpars(int &nx, double &xmi, double &xma, int &nq, double &qmi, double &qma, int &iord);
    void setlim( int ixmin, int iqmin, int iqmax, double dummy);
    void getlim( int iset, double &xmin, double &qmin, double &qmax, double &dummy);

    /************************************************************/
    /*  QCDNUM init routines from usrini.f                      */
    /************************************************************/

    void qcinit(int lun, string fname);
    void setlun(int lun, string fname);
    int nxtlun(int lmin);
    void setval(string opt, double val);
    void getval(string opt, double &val);
    void setint(string opt, int ival);
    void getint(string opt, int &ival);
    void qstore(string opt, int ival, double &val);

    /************************************************************/
    /*  QCDNUM parameter routines from usrparams.f              */
    /************************************************************/

    void setord(int iord);
    void getord(int &iord);
    void setalf(double as0, double r20);
    void getalf(double &as0, double &r20);
    void setcbt(int nfix, int iqc, int iqb, int iqt);
    void mixfns(int nfix, double r2c, double r2b, double r2t);
    void getcbt(int &nfix, double &q2c, double &q2b, double &q2t);
    void setabr(double ar, double br);
    void getabr(double &ar, double &br);
    int nfrmiq(int iset, int iq, int &ithresh);
    int nflavs(int iq, int &ithresh);
    void cpypar(double *array, int n, int iset);
    void usepar(int iset);
    int keypar(int iset);
    int keygrp(int iset, int igrp);
    void pushcp();
    void pullcp();

    /************************************************************/
    /*  QCDNUM pdf output routines from usrpdf.f                */
    /************************************************************/

    void allfxq(int iset, double x, double qmu2, double *pdfs, int n, int ichk);
    void allfij(int iset, int ix, int iq, double *pdfs, int n, int ichk);
    double bvalxq(int iset, int id, double x, double qmu2, int ichk);
    double bvalij(int iset, int id, int ix, int iq, int ichk);
    double fvalxq(int iset, int id, double x, double qmu2, int ichk);
    double fvalij(int iset, int id, int ix, int iq, int ichk);
    double sumfxq(int iset, double *c, int isel, double x, double qmu2, int ichk);
    double sumfij(int iset, double *c, int isel, int ix, int iq, int ichk);
    double fsplne(int iset, int id, double x, int iq);
    double splchk(int iset, int id, int iq);
  void fflist(int iset, double *c, int m, double *x, double *q, double *f, int n, int ichk);
  void fftabl(int iset, double *c, int isel, double *x, int nx, double *q, int nq, double *table, int m, int ichk);
  void ftable(int iset, double *c, int m, double *x, int nx, double *q, int nq, double *table, int ichk);
    void ffplot(string fnam, double(* fun)(int*, double*, bool*), int m, double zmi, double zma, int n, string txt);
    void fiplot(string fnam, double(* fun)(int*, double*, bool*), int m, double* zval, int n, string txt);

    /************************************************************/
    /*  QCDNUM workspace routines from usrstore.f               */
    /************************************************************/

    void maketab(double *w, int nw, int *jtypes, int npar,
                 int newt, int &jset, int &nwords);
    void setparw(double *w, int jset, double *par, int n);
    void getparw(double *w, int jset, double *par, int n);
    void dumptab(double *w, int jset, int lun,
                 string fnam, string fkey);
    void readtab(double *w, int nw, int lun, string fnam,
                 string fkey, int newt, int &jset,
                 int &nwords, int &ierr);

    /************************************************************/
    /*  QCDNUM weight routines from usrwgt.f                    */
    /************************************************************/

    void fillwt(int itype, int &idmin, int &idmax, int &nwds);
    void dmpwgt(int itype, int lun, string fname);
    void readwt(int lun, string fname, int &idmin, int &idmax, int &nwds, int &ierr);
    void wtfile(int itype, string fname);
    void nwused(int &nwtot, int &nwuse, int &ndummy);

    /************************************************************/
    /*  HQSTF routines                                          */
    /************************************************************/

    int ihqvers();
    void hqwords(int &ntotal, int &nused);
    void hqparms(double *qmass, double &a, double &b);
    double hqqfrmu(double qmu2);
    double hqmufrq(double Q2);
    void hswitch(int iset);
    void hqstfun(int istf, int icbt, double *def, double *x, double *Q2, double *f, int n, int ichk);
    void hqfillw(int istf, double *qmass, double aq, double bq, int &nused);
    void hqdumpw(int lun, string fname);
    void hqreadw(int lun, string fname, int &nused, int &ierr);

    /************************************************************/
    /*  ZMSTF routines                                          */
    /************************************************************/

    int izmvers();
    void zmwords(int &ntotal, int &nused);
    void zmdefq2(double a, double b);
    void zmabval(double &a, double &b);
    double zmqfrmu(double qmu2);
    double zmufrmq(double Q2);
    void zswitch(int iset);
    void zmstfun(int istf, double *def, double *x, double *Q2, double *f, int n, int ichk);
    void zmfillw(int &nused);
    void zmdumpw(int lun, string fname);
    void zmreadw(int lun, string fname, int &nused, int &ierr);
    void zmwfile(string fname);
}

#endif
