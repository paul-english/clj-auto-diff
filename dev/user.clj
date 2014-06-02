(ns user
  (:use midje.repl)
  (:require [clojure.tools.namespace.repl :refer [refresh refresh-all]]))

(defn start []
  (autotest))

(defn stop []
  (autotest :pause))

(defn reset []
  (stop)
  (refresh :after 'user/start))
