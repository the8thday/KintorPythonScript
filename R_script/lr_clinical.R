#! /urs/bin/env Rscript
# Title     : logistic analysis for clinical research
# Objective : 临床预测模型，基于回归模型的影响因素分析
# Created by: congliu
# Created on: 2021/8/13

library(tidymodels)
library(easystats)

# data
data(biopsy, package = 'MASS')
biopsy.omit <- biopsy %>%
  drop_na() %>%
  select(-ID)

data(Sacramento, package = 'modeldata')
set.seed(42)
s_split <- initial_split(Sacramento, prop = 0.8)
X_train <- training(s_split)
X_test <- testing(s_split)

# 应用条件 LINE=线性（主要对于定量数据、对于分类数据自变量没有要求）、独立、正态分布方差qixing、方差齐性
# 对于独立性、正态分布，主要要求模型的残差符合要求，对于不满足条件的数据或对Y值进行转换或用线性混合模型

# model
lm_model <- parsnip::linear_reg(
  mode = 'regression',
  penalty = NULL
) %>%
  set_engine('lm')

lm_recipe <- recipe(price ~ beds + baths + sqft + city, data = X_train) %>%
  step_dummy(city)

lm_wf <- workflow() %>%
  add_model(lm_model) %>%
  add_formula(price ~ beds + baths + sqft)

lm_fit <- lm_model %>% fit(price ~ beds + baths + sqft,
                 data = X_train)

lm_last_fit <- lm_wf %>% last_fit(
  split = s_split
)

# 模型条件诊断 和 总体拟合效果评价
lm_fit %>% performance::check_model()
lm_fit %>% performance::check_outliers()
r2(lm_fit)

lm_fit %>% model_performance()

names(lm_fit) #parsnip object
glance(lm_fit) # 查看模型performance
broom::tidy(lm_fit)
residuals(lm_fit$fit)
vip::vip(lm_fit) # 特征重要性 based F-statistic

# 预测新数据
X_test_res <- predict(lm_fit, new_data = X_test) %>%
  bind_cols(X_test)

yardstick::rmse(X_test_res,
     truth = price,
     estimate = `.pred`
)
yardstick::rsq(X_test_res,
     truth = price,
     estimate = `.pred`
)
# although we can use metrics set
multi_metric <- metric_set(rmse, rsq, mae)
multi_metric(X_test_res, truth = price, estimate = `.pred`)


# R2 plot
ggplot(data = X_test_res,
       mapping = aes(x = .pred, y = price)) +
  geom_point(color = '#006EA1') +
  geom_abline(intercept = 0, slope = 1, color = 'orange') +
  labs(title = 'Linear Regression Results',
       x = 'Predicted price',
       y = 'Actual price') +
  theme_bw()


# 对于 last_fit 的performance 操作
lm_last_fit %>% collect_metrics()

lm_recipe %>%
  prep() %>%
  bake(new_data = X_test) # if we got the correct transformations?

lm_last_res <- lm_last_fit %>% collect_predictions()
