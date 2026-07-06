# Протокол верификации списка литературы (2026-07-05)

Каждая позиция проверена веб-поиском по первоисточникам (IEEE Xplore, MDPI,
Google Patents, DBLP, arXiv): название, авторы, том, номер, страницы, год,
DOI. Полные журналы агентов — в workflow `verify-references`.

| № | Статус | Что исправлено / подтверждено |
|---|---|---|
| 1 | CORRECT | Mahony et al., IEEE TAC 53(5):1203–1218, 2008; DOI 10.1109/TAC.2008.923738 |
| 2 | CORRECTED (оформление) | ICORR 2011, Zurich, pp. 1–7; DOI 10.1109/ICORR.2011.5975346; не путать с неопубликованным отчётом 2010 г. |
| 3 | CORRECT | Valenti et al., Sensors 15(8):19302–19330, 2015; DOI 10.3390/s150819302 |
| 4 | CORRECT | VQF, Information Fusion 91:187–204, 2023; DOI 10.1016/j.inffus.2022.10.014 |
| 5 | CORRECT | BROAD, Data 6(7):72, 2021; DOI 10.3390/data6070072 (в названии em-dash) |
| 6 | CORRECT | Sabatini, IEEE TBME 53(7):1346–1356, 2006; DOI 10.1109/TBME.2006.875664 |
| 7 | **CORRECTED** | Верный порядок авторов: **Afzal, Renaudin, Lachapelle**; название "Use of Earth's magnetic field for mitigating gyroscope errors regardless of magnetic perturbation", Sensors 11(12):11390–11414, 2011; DOI 10.3390/s111211390. Важно: термин «MARU» введён в [8] (MAGYQ, 2014); [7] — техника квазистатического поля (Ḃ=−ω×B). Текст статьи скорректирован соответственно. «Complete Triaxis Magnetometer Calibration…» — ДРУГАЯ работа (J. Sensors 2010), не цитируется |
| 8 | CORRECT | MAGYQ, Sensors 14(12):22864–22890, 2014; DOI 10.3390/s141222864; вводит термины MARU и AGU |
| 9 | CORRECT | Roetenberg et al., IEEE TNSRE 13(3):395–405, 2005; DOI 10.1109/TNSRE.2005.847353 |
| 10 | CORRECT | Suh, IEEE TIM 59(12):3296–3305, 2010; DOI 10.1109/TIM.2010.2047157 |
| 11 | **CORRECTED** | Журнальная версия найдена: Fang, Haile, Wang, "Robust extended Kalman filtering for systems with measurement outliers", IEEE TCST 30(2):795–802, 2022; DOI 10.1109/TCST.2021.3077535. Конференц-версия: CDC 2018, pp. 6390–6395; DOI 10.1109/CDC.2018.8619140. VISKF: Sun et al., Micromachines 16(9):1036, 2025; DOI 10.3390/mi16091036 |
| 12 | CORRECTED (полное название) | Foxlin, US 6 361 507 B1, 26.03.2002; полное название "...for tracking human head and other similarly sized body"; правообладатель MIT |
| 13 | CORRECTED (добавлен изобретатель) | M. D. Fortier, US 8 645 063 B2, 04.02.2014; Custom Sensors & Technologies |
| 14 | CORRECT | Kok & Schön, IEEE SPL 26(11):1673–1677, 2019; DOI 10.1109/LSP.2019.2943995 |
| 15 | CORRECTED (подтверждена уместность) | Hua, Ducard, Hamel, Mahony, Rudin, IEEE TCST 22(1):201–213, 2014; DOI 10.1109/TCST.2013.2251635. Аннотация явно заявляет «measurement decoupling strategy... roll and pitch estimation robust to magnetic disturbances» — ссылка уместна для клейма о развязке каналов |
| 16 | CORRECT | github.com/xioTechnologies/Fusion, MIT-лицензия, активен; атрибуция S. Madgwick / x-io корректна (revised AHRS из гл. 7 диссертации Madgwick) |

## Карта цитирований (где какой [N] используется и что подтверждает)

| [N] | Разделы | Что подтверждает |
|---|---|---|
| [1] | Введение, 2.1, (2.4 неявно) | фильтр Mahony: PI-связь на SO(3), интегральный bias |
| [2] | Введение, 2.1, 2.4, 4 | фильтр Madgwick: градиентный шаг постоянной величины β |
| [3] | Введение (п.3), 2.1 | Valenti: развязка tilt/heading, адаптивный вес по норме |
| [4] | Введение, 2.2–2.5, 3, 6.4, 7 | VQF: архитектура, таймерный автомат, SOTA на BROAD, Cortex-M4 |
| [5] | 2.3, 5, 6.2 | датасет BROAD, официальные метрики и протокол TAGP |
| [6] | Введение, 2.2, 2.5 | Sabatini: кватернионный EKF, пороговые тесты поля |
| [7] | 2.5 | техника квазистатического поля Ḃ=−ω×B (Afzal et al.) |
| [8] | 2.5 | MAGYQ: термины MARU/AGU, инновации из теоремы переноса |
| [9] | 2.5 | Roetenberg: компенсация магнитных помех, онлайн-модель помехи |
| [10] | Введение, 2.2 | Suh: адаптивный косвенный КФ, инфляция ковариаций |
| [11] | Введение, 2.2, 4 | насыщение инновации: bounded error при bounded outliers |
| [12] | 2.4 | Foxlin: slew-rate-ограниченная компенсация дрейфа (прототип) |
| [13] | 2.4 | статические rate-лимиты коррекций (прототип, действующий патент) |
| [14] | 2.4 | Kok & Schön: градиентный шаг фиксированной длины |
| [15] | Введение (п.3) | развязка roll/pitch от магнитометра в нелинейном наблюдателе |
| [16] | Введение (п.1), 2.3, 7, 8 | Fusion: гейтирование + recovery-счётчики + bias в статике |

Патент US 10 852 846 упоминается в тексте (разд. 2.5) без номера в списке —
как патентный риск; при подаче можно добавить отдельной позицией.
