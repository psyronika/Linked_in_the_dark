# NL-Telegram-nets

This repository contains, data cleaning, processing, BERTopic modelling, network analysis and network visualisations of Telegram chats/channels that were used for the manuscript:

### Simon, M., Welbers, K., Kroon, A, & Trilling, D. Linked in the Dark: A network approach to understanding information flows within the Dutch Telegramsphere (forthcoming).

#### Abstract

Recent studies have shown that the stricter content moderation policies imposed by mainstream social networking sites (SNSs) stimulated the growth of low-moderated but relatively open discussion platforms such as Telegram. Despite Telegram’s growing popularity among (deplatformed) digital exiles, and high potential for news dissemination, information consumption, mobilization, and radicalization, little is known about information flows with respect to politically and socially relevant topics within the Telegramsphere. We scrutinize the Telegramsphere as an information-sharing ecosystem of current affairs by uncovering how information flows indicated by content-overlap and shared users influenced the structure of Telegram networks and shaped communities over time. Using state-of-the-art web-mining, neural topic modelling, and social network analysis techniques on a unique data set that spans the full messaging history (N = 2, 033, 661) of 174 Dutch-language public Telegram chats/channels, we show that over time, conspiracy-themed, far-right activist, and COVID-19-sceptical communities dominated the Dutch Telegramsphere of current affairs. Our findings raise concerns with respect to Telegram’s polarisation and radicalization capacity in the context of consuming socially and politically relevant information online.

## Data collection

The data used for the analyses was collected via a python script that is available via https://github.com/psyronika/Telegram-scraper.

## In this repo

- R scripts to conduct [network analysis] and network visualizations (https://github.com/psyronika/NL-Telegram-nets/blob/main/src/analysis)
- Python notebooks to conducet BERTopic modelling:
  - [prepare data for BERTopic](https://github.com/psyronika/NL-Telegram-nets/blob/main/src/data-processing/Prepare_data_for_BERTopic.ipynb)
  - [fine-tune BERTopic model](https://github.com/psyronika/NL-Telegram-nets/tree/main/src/analysis#:~:text=Fine_tune_BERTopic.ipynb)
  - [topics per class](https://github.com/psyronika/NL-Telegram-nets/blob/main/src/analysis/Topics_per_class.ipynb)
  - [topics over time](https://github.com/psyronika/NL-Telegram-nets/blob/main/src/analysis/Topics_over_time.ipynb)
- [Figures from the paper](https://github.com/psyronika/NL-Telegram-nets/tree/main/figures)



## More information
More information about the team and the project can be found via https://newsflows.eu/ 

#### This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (Grant agreement No. 947695). 
