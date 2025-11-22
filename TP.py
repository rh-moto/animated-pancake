# Load modules
import numpy as np
import pandas as pd
import os
from scipy.optimize import minimize
from scipy import integrate
from datetime import datetime
from xbbg import blp
from datetime import timedelta
from tqdm import tqdm

# Define functions
def daily_data_update_JGB(start_date=datetime(2000,1,1),end_date =None):
    JGB_ticker = [
                "JGBS1 Index",
                  "JGBS2 Index",
                  "JGBS3 Index",
                  "JGBS4 Index",
                  "JGBS5 Index",
                  "JGBS6 Index",
                  "JGBS7 Index",
                  "JGBS8 Index",
                  "JGBS9 Index",
                  "JGBS10 Index",
                  "JGBS15 Index",
                  "JGBS20 Index",
                  "JGBS25 Index",
                  "JGBS30 Index",
                  "JGBS40 Index"]
    #start_date = datetime(2024,1,1)
    start_y, start_m, start_d = start_date.year,start_date.month,start_date.day
    blp_start_date = str(start_y)+"-"+str(start_m)+"-"+str(start_d)
    if end_date == None:
      end_date = datetime.today() - timedelta(1)
      end_y, end_m, end_d = end_date.year, end_date.month, end_date.day
    blp_end_date = str(end_y)+"-"+str(end_m)+"-"+str(end_d)
    for  i  in range(len(JGB_ticker)):
      ticker = JGB_ticker[i]
      tmp = blp.bdh(tickers=ticker, flds=['last_price'],start_date= blp_start_date, end_date=blp_end_date)
      tmp = tmp[ticker].reset_index()
      if ticker == "MUTKCALM Index":
        tmp.columns = ["date","TONA"]
      else:
        tmp.columns = ["date",ticker]
      if i ==0:
        res  = tmp
      else:
        res  = pd.merge(res,tmp,on="date",how="outer" )
    return res

def daily_data_update_forwadswap_JPY(fwd_tenor,start_date=datetime(2000,1,1),end_date =None):
    fwd_end = ["1M","3M","6M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","20Y","30Y"]
    JGB_ticker = [f"S0195FS {fwd_tenor}{y} BLC Curncy" for y in fwd_end]
    start_y, start_m, start_d = start_date.year,start_date.month,start_date.day
    blp_start_date = str(start_y)+"-"+str(start_m)+"-"+str(start_d)
    if end_date == None:
      end_date = datetime.today() - timedelta(1)
      end_y, end_m, end_d = end_date.year, end_date.month, end_date.day
    blp_end_date = str(end_y)+"-"+str(end_m)+"-"+str(end_d)
    for  i  in range(len(JGB_ticker)):
      ticker = JGB_ticker[i]
      tmp = blp.bdh(tickers=ticker, flds=['last_price'],start_date= blp_start_date, end_date=blp_end_date)
      tmp = tmp[ticker].reset_index()
      if ticker == "MUTKCALM Index":
        tmp.columns = ["date","TONA"]
      else:
        tmp.columns = ["date",ticker]
      if i ==0:
        res  = tmp
      else:
        res  = pd.merge(res,tmp,on="date",how="outer" )
    return res

def daily_data_update_for_term_premium_est(start_date=datetime(2000,1,1),end_date =None):
    fwd_tenors =  ["1M","3M","6M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","15Y","20Y","30Y"]
    for fwd in fwd_tenors:
      tmp = daily_data_update_forwadswap_JPY(fwd,start_date=datetime(2000,1,1),end_date =None)
      dir_path =  "Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\input"
      os.chdir(dir_path)
      tmp.to_excel(f"ForwardSwap{fwd}.xlsx",index=False)
      print(f"forward swap data updated ({fwd})")
    JGB = daily_data_update_JGB(start_date=datetime(2000,1,1),end_date =None)
    JGB.to_excel("JGB_spot.xlsx",index=False)
    print("JGB data updated")
    return

def nelson_siegel(tau, beta0, beta1, beta2, lambda_):
    term1 = beta0
    term2 = beta1 * (1 - np.exp(-tau / lambda_)) / (tau / lambda_)
    term3 = beta2 * ((1 - np.exp(-tau / lambda_)) / (tau / lambda_) - np.exp(-tau / lambda_))
    return term1 + term2 + term3

def coef_nelson_siegel(tau,lambda_):
    term1 = 1
    term2 = (1 - np.exp(-tau / lambda_)) / (tau / lambda_)
    term3 =  ((1 - np.exp(-tau / lambda_)) / (tau / lambda_) - np.exp(-tau / lambda_))
    return np.array([term1,term2,term3])


# 誤差関数（最小化する対象）
def error_function(params, tau, yields,lambda_=8 ):
    beta0, beta1, beta2 = params
    model_yields = nelson_siegel(tau, beta0, beta1, beta2, lambda_)
    return np.sum((yields - model_yields)**2)

# 決定係数 R^2 の計算
def calculate_r_squared(y_obs, y_pred):
    ss_res = np.sum((y_obs - y_pred) ** 2)
    ss_tot = np.sum((y_obs - np.mean(y_obs)) ** 2)
    return 1 - (ss_res / ss_tot)

# フィッティング関数 (引数は yield_data のみ)
def fit_nelson_siegel(yield_data, maturities):
    results = []  # 各日のパラメータを保存するリスト
    lambda_value = 8
    # 各日について最適化を実行
    for i in range(yield_data.shape[0]):
        yields = yield_data.iloc[i].values  # pandas DataFrameから利回りデータを取得
        # 初期推定値
        initial_params = [0.02, -0.02, 0.01]  # [beta0, beta1, beta2, lambda]
        # 制約付き最適化
        res = minimize(error_function, initial_params, args=(maturities, yields,lambda_value),
                       method='SLSQP', bounds=[(-np.inf, np.inf), (-np.inf, np.inf), (-np.inf, np.inf)])
        if res.success:
            beta0, beta1, beta2 = res.x
            lambda_=lambda_value
            # モデルによる予測利回り
            model_yields = nelson_siegel(maturities, beta0, beta1, beta2, lambda_)
            # 決定係数 R^2 を計算
            r_squared = calculate_r_squared(yields, model_yields)
            results.append([beta0, beta1, beta2, lambda_, r_squared])
        else:
            results.append([np.nan] * 5)  # 失敗時はNaNを入れる
    # 結果をデータフレームに変換
    return pd.DataFrame(results, columns=['Beta0', 'Beta1', 'Beta2', 'Lambda', 'R_squared'])


def NS_data_fwd():
    dir_path =   "Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\input"
    os.chdir(dir_path)
    dic_tenor = dict({"1M":1/12,"3M":3/12,"6M":6/12,"1Y":1,"2Y":2,"3Y":3,"4Y":4,"5Y":5,"6Y":6,"7Y":7,"8Y":8,"9Y":9,"10Y":10,"15Y":15,"20Y":20,"30Y":30})
    fwd_tenors = ["1M","3M","6M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","20Y","30Y"]
    filenames = [f"ForwardSwap{x}.xlsx" for x in fwd_tenors]
    NS_res = []
    for j,x in enumerate(fwd_tenors):
        df1 = pd.read_excel(filenames[j],header=0)
        temp_df = pd.DataFrame({"date":df1.iloc[:,0]})
        for c in df1.columns:
            if "date"  not in c:
                temp_df[c] = df1[c]
        #temp_df.to_excel(f"NS_params_FWD_{x}.xlsx",index=False)
        NS_res.append(temp_df)
    print(f"FWD swap data {fwd_tenors} years ahead")
    return NS_res

def NS_para_fwd():
    NS_res = NS_data_fwd()
    maturities = np.array([1/12,3/12,6/12,1,2,3,4,5,6,7,8,9,10,20,30])
    fwd_tenors = ["1M","3M","6M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","20Y","30Y"]
    fwd_end = ["1M","3M","6M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","20Y","30Y"]
    NS_param_res =[]
    for i,x in enumerate(NS_res):
        tmp = fit_nelson_siegel(x.iloc[:,1:], maturities)
        tmp["date"] = x["date"]
        tmp.to_excel(f"Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\input\\NS_params\\NS_params_FWD_{fwd_tenors[i]}.xlsx",index=False)
        NS_param_res.append(tmp)
    print(f"Estimated the NS parameters in the FWD swap data {fwd_tenors} years ahead {fwd_end} years")
    return NS_param_res

def load_NS_para(tgt_dt):
    res= pd.DataFrame({})
    fwd_tenors = ["1M","3M","6M","1Y","2Y","3Y"]#"6Y","8Y","9Y",,"4Y","5Y","7Y","10Y","20Y","30Y"
    maturities = [1/12,3/12,6/12,1,2,3]
    for i,x in enumerate(fwd_tenors):
        tmp = pd.read_excel(f"Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\input\\NS_params\\NS_params_FWD_{fwd_tenors[i]}.xlsx",header=0)
        res = pd.concat([res,tmp.loc[tmp["date"] == tgt_dt]],axis=0)
    res["tenor"] = maturities
    res= res.reset_index(drop=True)
    res["short_rate_path"] = nelson_siegel(1/365, res["Beta0"], res["Beta1"], res["Beta2"], res["Lambda"])
    return res

def get_NS_params_short_rate_P_measure(data):
    res = fit_nelson_siegel(pd.DataFrame({"y":list(data["short_rate_path"])}).T, np.array(list(data["tenor"])))
    return list(res.values[0])

def get_NS_params_integrated_vol(data):
    integ_vol = ((data["short_rate_path"].diff().dropna() )**2).cumsum()
    integ_vol_tenor = np.array(list(data["tenor"])[1:])
    res = fit_nelson_siegel(pd.DataFrame({"y":list(integ_vol)}).T, integ_vol_tenor)
    return list(res.values[0])

def get_fwd_spot_rate(tgt_date,tenor):
    df_tmp = load_NS_para(tgt_date)[0]
    #tenor = np.array([1,2,5,7,10])#
    df_spot = pd.read_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\input\JGB_spot.xlsx",header=0)
    spot_rate = df_spot.loc[df_spot["date"]==tgt_date][[f"JGBS{_} Index" for _ in tenor]].reset_index(drop=True)
    spot_rate = np.array(spot_rate)[0]
    future_policy_rate = np.array(df_tmp.loc[(df_tmp["tenor"].apply(lambda x: x in list(tenor)) )]["short_rate_path"])
    return [spot_rate,future_policy_rate]


def get_NS_params_integrated_vol(data):
    integ_vol = ((data["short_rate_path"].diff().dropna() )**2).cumsum()
    integ_vol_tenor = np.array(list(data["tenor"])[1:])
    res = fit_nelson_siegel(pd.DataFrame({"y":list(integ_vol)}).T, integ_vol_tenor)
    return list(res.values[0])


def est_short_rate_risk_neutral(sr_params,vol_params):
    def er(params):
        y = np.array([ nelson_siegel(i*0.1,sr_params[0],sr_params[1],sr_params[2],sr_params[3]) for i in range(1,301)])
        vol = np.array([ nelson_siegel(i*0.1,vol_params[0],vol_params[1],vol_params[2],vol_params[3]) for i in range(1,301)])
        y_Q = np.array([ nelson_siegel(i*0.1,params[1],params[2],params[3],8) for i in range(1,301)])
        diff = y -vol*params[0] -y_Q
        return sum(diff**2)
    initial_params = [0.1,0.1,0.1,0.1]
    optim_res = minimize(er, initial_params, #args=(spot_rate, future_policy_rate,tenor,LAM),
                       method='SLSQP',bounds=[(-np.inf, np.inf), (-np.inf, np.inf), (-np.inf, np.inf), (-np.inf, np.inf)])#constraints=linear_constraint
    return optim_res.x

def estimate_short_rate_path_risk_neutral(tgt_date):
    data=  load_NS_para(tgt_date)
    vol_params = get_NS_params_integrated_vol(data)
    sr_params = get_NS_params_short_rate_P_measure(data)
    rp_para = est_short_rate_risk_neutral(sr_params,vol_params)
    #sr = np.array([nelson_siegel(i,rp_para[1],rp_para[2],rp_para[3],8) for i in range(1,30)])
    res= pd.DataFrame({"date":[tgt_date]})
    res["coef_market_price_risk"] = rp_para[0]
    res["beta0_risk_neutral_measure"] = rp_para[1]
    res["beta1_risk_neutral_measure"] = rp_para[2]
    res["beta2_risk_neutral_measure"] = rp_para[3]
    res["lambda_risk_neutral_measure"] = 8
    res["beta0_real_measure"] = sr_params[0]
    res["beta1_real_measure"] = sr_params[1]
    res["beta2_real_measure"] = sr_params[2]
    res["lambda_risk_real_measure"] = 8
    res["beta0_integrated_vol"] = vol_params[0]
    res["beta1_integrated_vol"] = vol_params[1]
    res["beta2_integrated_vol"] = vol_params[2]
    res["lambda_integrated_vol"] = 8
    for t in [2,5,10,20,30]:
        res[f"RNR_{t}y"] = (1/t)*integrate.quad(nelson_siegel,0,t,args=(rp_para[1], rp_para[2], rp_para[3], 8))[0]
    return res

def estimate_RNR():
    tmp = pd.read_excel(f"Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\input\\NS_params\\NS_params_FWD_1M.xlsx",header=0)
    tmp =tmp.loc[tmp["date"]>=datetime(2010,1,1)].reset_index(drop=True)
    date_list = list(tmp["date"])
    if os.path.exists("Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\output\\params_future_policy_rate_curve.xlsx"):
        old_data = pd.read_excel("Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\output\\params_future_policy_rate_curve.xlsx",header=0)
        old_data["date"] = pd.to_datetime(old_data["date"] )
        date_list = [d for d in date_list if d > max(old_data["date"] )]
    else:
        date_list = date_list
    res = pd.DataFrame({})
    for d in tqdm(date_list):
        res = pd.concat([res,estimate_short_rate_path_risk_neutral(d)],axis =0)
    if os.path.exists("Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\output\\params_future_policy_rate_curve.xlsx"):
        res =pd.concat([old_data,res],axis=0).drop_duplicates(["date"]).sort_values(["date"])
    res.to_excel("Y:\\TKY\ops\\economics\\FXR\\Python\\Codes\\Term premium_revised\\output\\params_future_policy_rate_curve.xlsx",index=False)
    return


def get_term_premium_estimation():
    df_jgb = pd.read_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\input\JGB_spot.xlsx",header =0)
    df_rnr = pd.read_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\output\params_future_policy_rate_curve.xlsx",header =0)
    df_jgb =df_jgb[["date","JGBS2 Index","JGBS5 Index","JGBS10 Index","JGBS20 Index","JGBS30 Index"]]
    df_jgb =df_jgb.loc[df_jgb["date"]>=datetime(2010,1,1)].reset_index(drop=True)
    df_jgb.columns = ["date","JGB_2y","JGB_5y","JGB_10y","JGB_20y","JGB_30y"]
    df_rnr =df_rnr.loc[df_rnr["date"]>=datetime(2010,1,1)].reset_index(drop=True)
    df_rnr = df_rnr[["date","RNR_2y","RNR_5y","RNR_10y","RNR_20y","RNR_30y"]]
    df = pd.merge(df_jgb,df_rnr,on="date",how="left")
    for i in [2,5,10,20,30]:
        df[f"TP_{i}y"] = df[f"JGB_{i}y"] -df[f"RNR_{i}y"]
    df = df.dropna().reset_index(drop=True)
    df.to_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\output\JGB_term_premium_estimation.xlsx",index =0)
    df["date"] = pd.to_datetime(df["date"])
    res_df = df.loc[(df["date"]<=datetime(2012,5,29))|(df["date"]>=datetime(2012,6,8))]
    res_df.to_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\output\JGB_term_premium_estimation_outlier_truncated.xlsx",index =0)
    return

def add_40y_term_premium():
    df_jgb = pd.read_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\input\JGB_spot.xlsx",header =0)
    df_tp = pd.read_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\output\JGB_term_premium_estimation.xlsx",header =0)
    df_jgb = df_jgb[["date","JGBS40 Index"]]
    df_jgb.columns = ["date","JGB_40y"]
    df_tp = pd.merge(df_tp,df_jgb,on="date",how="left")
    df_tp["RNR_40y"] = (30*df_tp["RNR_30y"] -20*df_tp["RNR_20y"] +30*df_tp["RNR_30y"])/40
    df_tp["TP_40y"] = df_tp["JGB_40y"] - df_tp["RNR_40y"]
    df_tp= df_tp[["date","JGB_2y","JGB_5y","JGB_10y","JGB_20y","JGB_30y","JGB_40y","RNR_2y","RNR_5y","RNR_10y","RNR_20y","RNR_30y","RNR_40y","TP_2y","TP_5y","TP_10y","TP_20y","TP_30y","TP_40y"]]
    df_tp.to_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\output\JGB_term_premium_estimation.xlsx",index =0)
    df_tp = df_tp.loc[(df_tp["date"]<=datetime(2012,5,29))|(df_tp["date"]>=datetime(2012,6,8))]
    df_tp.to_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\output\JGB_term_premium_estimation_outlier_truncated.xlsx",index =0)
    return


def daily_data_update_TONA(start_date=datetime(2000,1,1),end_date =None):
    JGB_ticker = ["MUTKCALM Index"]
    start_y, start_m, start_d = start_date.year,start_date.month,start_date.day
    blp_start_date = str(start_y)+"-"+str(start_m)+"-"+str(start_d)
    if end_date == None:
      end_date = datetime.today() - timedelta(1)
      end_y, end_m, end_d = end_date.year, end_date.month, end_date.day
    blp_end_date = str(end_y)+"-"+str(end_m)+"-"+str(end_d)
    for  i  in range(len(JGB_ticker)):
      ticker = JGB_ticker[i]
      tmp = blp.bdh(tickers=ticker, flds=['last_price'],start_date= blp_start_date, end_date=blp_end_date)
      tmp = tmp[ticker].reset_index()
      if ticker == "MUTKCALM Index":
        tmp.columns = ["date","TONA"]
      else:
        tmp.columns = ["date",ticker]
      if i ==0:
        res  = tmp
      else:
        res  = pd.merge(res,tmp,on="date",how="outer" )
    df_jgb = pd.read_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\input\JGB_spot.xlsx",header =0)
    df_jgb = df_jgb[["date"]]
    res["date"] = pd.to_datetime(res["date"])
    res = pd.merge_asof(df_jgb,res)
    res =res.loc[res["date"]>=datetime(2010,1,1)].reset_index(drop=True)
    res.to_excel(r"\\intranet.barcapint.com\dfs-apac\Group\TKY\ops\economics\FXR\Python\Codes\Term premium_revised\input\TONA.xlsx",index =0)
    return res

# run estimation
daily_data_update_JGB()
daily_data_update_for_term_premium_est()
NS_para_fwd()
estimate_RNR()
get_term_premium_estimation()
daily_data_update_TONA()
add_40y_term_premium()