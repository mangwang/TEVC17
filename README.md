### Source Code for Rerunning all the Experiments of my TEVC Paper ###


1. The **_main.m_** file is the entry point to run all the experements;

2. The first experiment is visulization of generational measures with convergence curve. In this experiment, the parameters _limFEs_ and _M_ are separately set to _maxFEs/5_ and _5_;

3. The second experiment performs the algorithm selection task via the proposed framework based on the _evp(P_i)_ sequences. Note that the _uniFEs_ should be calculated before estimating the measures, while all the measures are already estimated in the first experiment, here we just reuse them again and unify the budget to _uniFEs_;

4. The statistical results of the proposed approach and performance evaluation approach are saved in _result_ folder, separately named _evp.csv_ and _bench.csv_. Note that the statistically significant data is shown followed by (*);

5. The benchmark results are pre-excuted and saved in the _result/benchmark_ folder.

If any questions, please no hesitate to email me (mangwang AT mail.ustc.edu.cn).