(ns ad.core
  (:refer-clojure :exclude [* + - / == < > <= >= zero? number? pos? neg?])
  (:require [clojure.core.match :refer [match]]
            [clojure.core.matrix :as matrix]
            [clojure.core.matrix.operators :as operators :refer [/]]))

;; There should be a way to do this without the mutable state here.
;; This is just how it was done in the scheme lib.
(def e (atom 0))

(def <_e clojure.core/<)

(deftype DualNumber [epsilon primal perturbation])

(definterface ITape
  (getFanout [])
  (setFanout [val])
  (getSensitivity [])
  (setSensitivity [val]))

(deftype Tape [epsilon primal factors tapes
               ^:volatile-mutable fanout
               ^:volatile-mutable sensitivity]
  ITape

  (getFanout [_] fanout)
  (setFanout [_ val] (set! fanout val))

  (getSensitivity [_] sensitivity)
  (setSensitivity [_ val] (set! sensitivity val)))

(defn dual-number? [n] (instance? DualNumber n))
(defn tape? [t] (instance? Tape t))

(defn ad-type [x]
  (cond (dual-number? x) :dual
        (tape? x) :tape))

(defn compare-differentiable [x1 x2]
  (try (cond (<_e (.epsilon x1) (.epsilon x2)) :<
             (<_e (.epsilon x2) (.epsilon x1)) :>
             :else nil)
       (catch Exception e nil)))

(defn new-tape [epsilon primal factors tapes]
  (Tape. epsilon primal factors tapes 0 0))

(defn tapify [x]
  (new-tape @e x '[] '[]))

(declare ^:dynamic **)
(declare * + cos cosh)

(defn lift-real->real [f df-dx]
  (letfn [(self [x]
            (match (ad-type x)
                   :dual (DualNumber. (.epsilon x) (self (.primal x)) (* (df-dx (.primal x)) (.perturbation x)))
                   :tape (new-tape (.epsilon x) (self (.primal x)) [(df-dx (.primal x))] [x])
                   _ (f x)))]
    self))

(defn lift-real*real->real [f df-dx1 df-dx2]
  (letfn [(self [x1 x2]
            (match [(ad-type x1) (ad-type x2) (compare-differentiable x1 x2)]
                   [:dual :dual :<] (DualNumber. (.epsilon x2) (self x1 (.primal x2)) (* (df-dx2 x1 (.primal x2)) (.perturbation x2)))
                   [:dual :dual :>] (DualNumber. (.epsilon x1) (self (.primal x1) x2) (* (df-dx1 (.primal x1) x2) (.perturbation x1)))
                   [:dual :dual  _] (DualNumber. (.epsilon x1) (self (.primal x1) (.primal x2)) (+ (* (df-dx1 (.primal x1) (.primal x2)) (.perturbation x1)) (* (df-dx2 (.primal x1) (.primal x2)) (.perturbation x2))))
                   [:dual :tape :<] (new-tape (.epsilon x2) (self x1 (.primal x2)) [(df-dx2 x1 (.primal x2))] [x2])
                   [:dual :tape :>] (DualNumber. (.epsilon x1) (self (.primal x1) x2) (* (df-dx1 (.primal x1) x2) (.perturbation x1)))
                   [:dual     _  _] (DualNumber. (.epsilon x1) (self (.primal x1) x2) (* (df-dx1 (.primal x1) x2) (.perturbation x1)))
                   [:tape :dual :<] (DualNumber. (.epsilon x2) (self x1 (.primal x2)) (* (df-dx2 x1 (.primal x2)) (.perturbation x2)))
                   [:tape :dual :>] (new-tape (.epsilon x1) (self (.primal x1) x2) [(df-dx1 (.primal x1) x2)] [x1])
                   [:tape :tape :<] (new-tape (.epsilon x2) (self x1 (.primal x2)) [(df-dx2 x1 (.primal x2))] [x2])
                   [:tape :tape :>] (new-tape (.epsilon x1) (self (.primal x1) x2) [(df-dx1 (.primal x1) x2)] [x1])
                   [:tape :tape  _] (new-tape (.epsilon x1) (self (.primal x1) (.primal x2)) [(df-dx1 (.primal x1) (.primal x2)) (df-dx2 (.primal x1) (.primal x2))] [x1 x2])
                   [:tape     _  _] (new-tape (.epsilon x1) (self (.primal x1) x2) [(df-dx1 (.primal x1) x2)] [x1])
                   [    _ :dual  _] (DualNumber. (.epsilon x2) (self x1 (.primal x2)) (* (df-dx2 x1 (.primal x2)) (.perturbation x2)))
                   [    _ :tape  _] (new-tape (.epsilon x2) (self x1 (.primal x2)) [(df-dx2 x1 (.primal x2))] [x2])
                   [    _     _  _] (f x1 x2)))]
    self))

(defn lift-real-n->real [f df-dx1 df-dx2]
  (fn [& xs]
    (if (nil? xs)
      (f)
      (reduce (lift-real*real->real f df-dx1 df-dx2) xs))))

(defn lift-real-n+1->real [f df-dx df-dx1 df-dx2]
  (fn [& xs]
    (cond (nil? xs) (f)
          (nil? (rest xs)) ((lift-real->real f df-dx) (first xs))
          :else (reduce (lift-real*real->real f df-dx1 df-dx2) xs))))

(defn primal* [x]
  (match (ad-type x)
         :dual (primal* (.primal x))
         :tape (primal* (.primal x))
         _ x))

(defn lift-real-n->boolean [f]
  (fn [& xs] (apply f (map primal* xs))))

(def + (lift-real-n->real operators/+
                          (fn [x1 x2] 1)
                          (fn [x1 x2] 1)))

(def - (lift-real-n+1->real operators/-
                            (fn [x] -1)
                            (fn [x1 x2] 1)
                            (fn [x1 x2] -1)))

(def * (lift-real-n->real operators/*
                          (fn [x1 x2] x2)
                          (fn [x1 x2] x1)))

(def div (lift-real-n+1->real /
                              (fn [x] (- (div (** x 2))))
                              (fn [x1 x2] (div x2))
                              (fn [x1 x2] (- (div x1 (** x2 2))))))

(def sqrt (lift-real->real #(Math/sqrt %)
                           (fn [x] (div 1 (* 2 (sqrt x))))))

(def exp (lift-real->real #(Math/exp %)
                          (fn [x] (exp x))))

(def log (lift-real->real #(Math/log %)
                          (fn [x] (div x))))

(def ^:dynamic ** (lift-real*real->real operators/**
                                        (fn [x1 x2] (* x2 (** x1 (- x2 1))))
                                        (fn [x1 x2] (* (log x1) (** x1 x2)))))

(def mmul (lift-real*real->real matrix/mmul
                                (fn [x1 x2] x1)
                                (fn [x1 x2] x1)))

(def sin (lift-real->real #(Math/sin %)
                          (fn [x] (cos x))))

(def cos (lift-real->real #(Math/cos %)
                          (fn [x] (- (sin x)))))

(def tan (lift-real->real #(Math/tan %)
                          (fn [x] (+ 1 (** (tan x) 2)))))

(def asin (lift-real->real #(Math/asin %)
                           (fn [x] (div 1 (sqrt (- 1 (** x 2)))))))

(def acos (lift-real->real #(Math/acos %)
                           (fn [x] (- (div 1 (sqrt (- 1 (** x 2))))))))

(def atan (lift-real->real #(Math/atan %)
                           (fn [x] (div 1 (+ 1 (* x x))))))

(def sinh (lift-real->real #(Math/sinh %)
                           (fn [x] (cosh x))))

(def cosh (lift-real->real #(Math/cosh %)
                           (fn [x] (sinh x))))

(def tanh (lift-real->real #(Math/tanh %)
                           (fn [x] (- 1 (** (tanh x) 2)))))

(def == (lift-real-n->boolean operators/==))

(def < (lift-real-n->boolean clojure.core/<))

(def > (lift-real-n->boolean clojure.core/>))

(def <= (lift-real-n->boolean clojure.core/<=))

(def >= (lift-real-n->boolean clojure.core/>=))

(def zero? (lift-real-n->boolean clojure.core/zero?))

(def pos? (lift-real-n->boolean clojure.core/pos?))

(def neg? (lift-real-n->boolean clojure.core/neg?))

(def number? (lift-real-n->boolean clojure.core/number?))

(def dot (lift-real-n->real matrix/dot
                            (fn [x1 x2] x1)
                            (fn [x1 x2] x1)))

(def norm (lift-real->real #(-> (dot % %) sqrt)
                                (fn [x] x)))

(def distance (lift-real*real->real #(-> (- %1 %2) norm)
                                    (fn [x1 x2] x1)
                                    (fn [x1 x2] x1)))

(defn write-real [x]
  (match (ad-type x)
         :dual (do (write-real (.primal x)) x)
         :tape (do (write-real (.primal x)) x)
         _ (do (println x) x)))

(defn forward-mode [map-independent map-dependent f x x-perturbation]
  ;; Based on R6RS-ad, thus doesn't support tangent vector mode
  (swap! e inc)
  (let [y-forward (f (map-independent (fn [x x-perturbation] (DualNumber. @e x x-perturbation))
                                      x
                                      x-perturbation))]
    (swap! e dec)
    [(map-dependent (fn [y-forward]
                      (if (or (not (dual-number? y-forward))
                              (<_e (.epsilon y-forward) @e))
                        y-forward
                        (.primal y-forward)))
                    y-forward)
     (map-dependent (fn [y-forward]
                      (if (or (not (dual-number? y-forward))
                              (<_e (.epsilon y-forward) @e))
                        0
                        (.perturbation y-forward)))
                    y-forward)]))

(defn diff [f]
  (fn [x]
    (second (forward-mode (fn [f x x-perturbation] (f x x-perturbation))
                          (fn [f y-forward] (f y-forward))
                          f x 1))))

(defn directional-derivative-vector-F [f]
  (fn [x x-perturbation]
    (second
     (forward-mode (fn [f x x-perturbation] (map f x x-perturbation))
                   (fn [f y-forward] (map f y-forward))
                   f
                   x
                   x-perturbation))))

(defn gradient-vector-F [f]
  (fn [x]
    (map (fn [i]
           ((diff (fn [xi] (f (assoc x i xi))))
            (nth x i)))
         (range (count x)))))

(defn determine-fanout! [tape]
  (.setFanout tape (+ (.fanout tape) 1))
  (when (= (.fanout tape) 1)
    (map determine-fanout! (.tapes tape))))

(defn initialize-sensitivity! [tape]
  (.setSensitivity tape 0)
  (.setFanout tape (- (.fanout tape) 1))
  (when (zero? (.fanout tape))
    (map initialize-sensitivity! (.tapes tape))))

(defn reverse-phase! [sensitivity tape]
  (.setSensitivity tape (+ (.sensitivity tape) sensitivity))
  (.setFanout tape (- (.fanout tape) 1))
  (when (zero? (.fanout tape))
    (let [sensitivity (.sensitivity tape)]
      (map
       (fn [factor tape] (reverse-phase! (* sensitivity factor) tape))
       (.factors tape)
       (.tapes tape)))))

(defn reverse-mode [map-independent
                    map-dependent
                    for-each-dependent1!
                    for-each-dependent2!
                    f
                    x
                    y-sensitivities]
  ;; needs work: We don't support providing the y-sensitivies
  ;;             (potentially incrementally) after computing the
  ;;             primal in the forward phase.
  (swap! e inc)
  (let [x-reverse (map-independent tapify x)
        y-reverse (f x-reverse)
        x-sensitivities (map (fn [y-sensitivity]
                               (for-each-dependent1!
                                (fn [y-reverse]
                                  (when (and (tape? y-reverse)
                                             (not (<_e (.epsilon y-reverse) @e)))
                                    (determine-fanout! y-reverse)
                                    (initialize-sensitivity! y-reverse)))
                                y-reverse)
                               (for-each-dependent2!
                                (fn [y-reverse y-sensitivity]
                                  (when (and (tape? y-reverse)
                                             (not (<_e (.epsilon y-reverse) @e)))
                                    (determine-fanout! y-reverse)
                                    (reverse-phase! y-sensitivity y-reverse)))
                                y-reverse
                                y-sensitivity)
                               (map-independent #(.sensitivity %) x-reverse))
                             y-sensitivities)]
    (swap! e dec)
    [(map-dependent
      (fn [y-reverse]
        (if (or (not (tape? y-reverse)) (<_e (.epsilon y-reverse) @e))
          y-reverse
          (.primal y-reverse)))
      y-reverse)
     x-sensitivities]))

(defn derivative-R [f]
  (fn [x]
    (first (second (reverse-mode
                    (fn [f x] (f x))
                    (fn [f y-reverse] (f y-reverse))
                    (fn [f y-reverse] (f y-reverse))
                    (fn [f y-reverse y-sensitivity] (f y-reverse y-sensitivity))
                    f
                    x
                    [1])))))

(defn gradient-vector-R [f]
    (fn [x]
      (first (second (reverse-mode
                      (fn [f x] (map f x))
                      (fn [f y-reverse] (f y-reverse))
                      (fn [f y-reverse] (f y-reverse))
                      (fn [f y-reverse y-sensitivity] (f y-reverse y-sensitivity))
                      f
                      x
                      [1])))))

(defn f-gradient-vector-vector-R [f]
    (fn [x]
      (let [result (reverse-mode
                    (fn [f x] (map (fn [x] (map f x)) x))
                    (fn [f y-reverse] (f y-reverse))
                    (fn [f y-reverse] (f y-reverse))
                    (fn [f y-reverse y-sensitivity] (f y-reverse y-sensitivity))
                    f
                    x
                    [1])]
        [(first result) (first (second result))])))
