(ns ad.core)

;; TODO we'll want this to work with core.matrix somehow...
;; that could get complex fast....
;; dealing with the dualnumber & tape types will make
;; matrices tough to do.

;; There should be a way to do this without the mutable state here.
;; This is just how it was done in the scheme lib.
(def e (atom 0))

(def <_e <)

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

(defn new-tape [epsilon primal factors tapes]
  (Tape. epsilon primal factors tapes 0 0))

(defn tapify [x]
  (new-tape @e x '[] '[]))

(declare d* d+ dcos)

(defn lift-real->real [f df-dx]
  ;;(println "lift-real->real" "f:" f "df-dx:" df-dx)
  (letfn [(self [x]
            (cond
             (dual-number? x) (DualNumber. (.epsilon x)
                                           (self (.primal x))
                                           (d* (df-dx (.primal x))
                                               (.perturbation x)))
             (tape? x) (new-tape (.epsilon x)
                                 (self (.primal x))
                                 [(df-dx (.primal x))]
                                 [x])
             :else (f x)))]
    self))

(defn lift-real*real->real [f df-dx1 df-dx2]
  ;;(println "lift-real*real->real" "f:" f "df-dx1:" df-dx1 "df-dx2:" df-dx2)
  (letfn [(self [x1 x2]
            ;;(println "-- self" x1 x2)
            ;;(println "self x1" (dual-number? x1) (tape? x1))
            ;;(println "self x2" (dual-number? x2) (tape? x2))
            (cond
             (dual-number? x1) (cond
                                (dual-number? x2) (cond
                                                   (<_e (.epsilon x1) (.epsilon x2)) (DualNumber. (.epsilon x2)
                                                                                                  (self x1 (.primal x2))
                                                                                                  (d* (df-dx2 x1 (.primal x2))
                                                                                                      (.perturbation x2)))
                                                   (<_e (.epsilon x2) (.epsilon x1)) (DualNumber. (.epsilon x1)
                                                                                                  (self (.primal x1) x2)
                                                                                                  (d* (df-dx1 (.primal x1) x2)
                                                                                                      (.perturbation x1)))
                                                   :else (DualNumber. (.epsilon x1)
                                                                      (self (.primal x1) (.primal x2))
                                                                      (d+ (d* (df-dx1 (.primal x1)
                                                                                      (.primal x2))
                                                                              (.perturbation x1))
                                                                          (d* (df-dx2 (.primal x1)
                                                                                      (.primal x2))
                                                                              (.perturbation x2)))))
                                (tape? x2) (if (<_e (.epsilon x1) (.epsilon x2))
                                             (new-tape (.epsilon x2)
                                                       (self x1 (.primal x2))
                                                       [(df-dx2 x1 (.primal x2))]
                                                       [x2])
                                             (DualNumber. (.epsilon x1)
                                                          (self (.primal x1) x2)
                                                          (d* (df-dx1 (.primal x1) x2)
                                                              (.perturbation x1))))
                                :else (DualNumber. (.epsilon x1)
                                                   (self (.primal x1) x2)
                                                   (d* (df-dx1 (.primal x1) x2)
                                                       (.perturbation x1))))
             (tape? x1) (cond
                         (dual-number? x2) (if (<_e (.epsilon x1) (.epsilon x2))
                                             (DualNumber. (.epsilon x2)
                                                          (self x1 (.primal x2))
                                                          (d* (df-dx2 x1 (.primal x2))
                                                              (.perturbation x2)))
                                             (new-tape (.epsilon x1)
                                                       (self (.primal x1) x2)
                                                       [(df-dx1 (.primal x1) x2)]
                                                       [x1]))
                         (tape? x2) (cond
                                     (<_e (.epsilon x1) (.epsilon x2)) (new-tape (.epsilon x2)
                                                                                 (self x1 (.primal x2))
                                                                                 [(df-dx2 x1 (.primal x2))]
                                                                                 [x2])
                                     (<_e (.epsilon x2) (.epsilon x1)) (new-tape (.epsilon x1)
                                                                                 (self (.primal x1) x2)
                                                                                 [(df-dx1 (.primal x1) x2)]
                                                                                 [x1])
                                     :else (new-tape (.epsilon x1)
                                                     (self (.primal x1) (.primal x2))
                                                     [(df-dx1 (.primal x1) (.primal x2))
                                                      (df-dx2 (.primal x1) (.primal x2))]
                                                     [x1 x2]))
                         :else (new-tape (.epsilon x1)
                                         (self (.primal x1) x2)
                                         [(df-dx1 (.primal x1) x2)]
                                         [x1]))
             :else (cond
                    (dual-number? x2) (DualNumber. (.epsilon x2)
                                                   (self x1 (.primal x2))
                                                   (d* (df-dx2 x1 (.primal x2))
                                                       (.perturbation x2)))
                    (tape? x2) (new-tape (.epsilon x2)
                                         (self x1 (.primal x2))
                                         [(df-dx2 x1 (.primal x2))]
                                         [x2])
                    :else (f x1 x2))))]
    self))

(defn lift-real-n->real [f df-dx1 df-dx2]
  ;;(println "lift-real-n->real" "f:" f "df-dx1:" df-dx1 "df-dx2:" df-dx2)
  (fn [& xs]
    (if (nil? xs)
      (f)
      (reduce (lift-real*real->real f df-dx1 df-dx2) xs))))

(defn lift-real-n+1->real [f df-dx df-dx1 df-dx2]
  ;;(println "lift-real-n+1->real" "f:" f "df-dx:" df-dx "df-dx1:" df-dx1 "df-dx2:" df-dx2)
  (fn [& xs]
    (cond (nil? xs) (f)
          (nil? (rest xs)) ((lift-real->real f df-dx) (first xs))
          :else (reduce (lift-real*real->real f df-dx1 df-dx2) xs))))

(defn primal* [x]
  (cond (dual-number? x) (primal* (.primal x))
        (tape? x) (primal* (.primal x))
        :else x))

(defn lift-real-n->boolean [f]
  ;;(println "lift-real-n->boolean f:" f)
  (fn [& xs] (apply f (map primal* xs))))

(def d+ (lift-real-n->real + (fn [x1 x2] 1) (fn [x1 x2] 1)))

(def d- (lift-real-n+1->real - (fn [x] -1) (fn [x1 x2] 1) (fn [x1 x2] -1)))

(def d* (lift-real-n->real * (fn [x1 x2] x2) (fn [x1 x2] x1)))

(def ddiv (lift-real-n+1->real
         /
         (fn [x] (d- (ddiv (d* x x))))
         (fn [x1 x2] (ddiv x2))
         (fn [x1 x2] (d- (ddiv x1 (d* x2 x2))))))

(def dsqrt (lift-real->real #(Math/sqrt %) (fn [x] (ddiv (d* 2 (dsqrt x))))))

(def dexp (lift-real->real #(Math/exp %) (fn [x] (dexp x))))

(def dlog (lift-real->real #(Math/log %) (fn [x] (ddiv x))))

(def dexpt
  (lift-real*real->real #(Math/pow %1 %2)
                        (fn [x1 x2] (d* x2 (dexpt x1 (d- x2 1))))
                        (fn [x1 x2] (d* (dlog x1) (dexpt x1 x2)))))

(def dsin (lift-real->real #(Math/sin %) (fn [x] (dcos x))))

(def dcos (lift-real->real #(Math/cos %) (fn [x] (d- (dsin x)))))

(def datan (lift-real->real #(Math/atan %) (fn [x] (ddiv 1 (d+ 1 (d* x x))))))

(def d= (lift-real-n->boolean =))

(def d< (lift-real-n->boolean <))

(def d> (lift-real-n->boolean >))

(def d<= (lift-real-n->boolean <=))

(def d>= (lift-real-n->boolean >=))

(def dzero? (lift-real-n->boolean zero?))

(def dpositive? (lift-real-n->boolean pos?))

(def dnegative? (lift-real-n->boolean neg?))

(def dreal? (lift-real-n->boolean number?))

(defn write-real [x]
  (cond (dual-number? x) (do (write-real (.primal x)) x)
        (tape? x) (do (write-real (.primal x)) x)
        :else (do (println x) x)))

(defn forward-mode [map-independent map-dependent f x x-perturbation]
  ;; Based on R6RS-ad, thus doesn't support tangent vector mode
  ;;(println "forward-mode")
  ;;(println "-- map-independent" map-independent)
  ;;(println "-- map-dependent" map-dependent)
  ;;(println "-- f" f)
  ;;(println "-- x" x)
  ;;(println "-- x-perturbation" x-perturbation)
  (swap! e inc)
  (let [y-forward (f (map-independent (fn [x x-perturbation] (DualNumber. @e x x-perturbation))
                                      x
                                      x-perturbation))]
    ;;(println "y-forward" y-forward)
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

(defn derivative-F [f]
  (fn [x]
    ;;(println "deriv-f fn")
    (second (forward-mode (fn [f x x-perturbation]
                            ;;(println "map-independent" f x x-perturbation)
                            (f x x-perturbation))
                          (fn [f y-forward]
                            ;;(println "map-dependent" f y-forward)
                            (f y-forward))
                          f x 1))))

(defn directional-derivative-list-F [f]
  (fn [x x-perturbation]
    (second (forward-mode (fn [f x x-perturbation] (map f x x-perturbation))
                          (fn [f y-forward] (map f y-forward))
                          f
                          x
                          x-perturbation))))

(defn for-each-n [f n]
  (loop [i 0] (when (< i n) (f i) (recur (+ i 1)))))

;; TODO remove vector specific functions
;; TODO probably don't need, clojures map works on vectors just as
;; well as lists
#_(defn map-vector [f v & vs]
  (let [u (make-vector (count v))]
    (for-each-n
     (fn [i]
       (vector-set!
        u i (apply f (vector-ref v i) (map (fn [v] (vector-ref v i)) vs))))
     (count v))
    u))

(defn directional-derivative-vector-F [f]
  (fn [x x-perturbation]
    (second
     (forward-mode (fn [f x x-perturbation] (map f x x-perturbation))
                   (fn [f y-forward] (map f y-forward))
                   f
                   x
                   x-perturbation))))

;; TODO probably recreating something already in core
(defn replace-ith [x i xi]
  (if (zero? i)
    (cons xi (rest x))
    (cons (first x) (replace-ith (rest x) (- i 1) xi))))

;; TODO probably recreating something already in core
(defn map-n [f n]
  (loop [result []
         i 0]
    (if (= i n) (reverse result) (recur (cons (f i) result) (+ i 1)))))

(defn gradien-list-F [f]
  (fn [x]
    (map-n
     (fn [i]
       ((derivative-F (fn [xi] (f (replace-ith x i xi)))) (nth x i)))
     (count x))))

(defn replace-ith-vector [x i xi]
    (map (fn [j] (if (= j i) xi (nth x j)))
         (range (count x))))

(defn gradient-vector-F [f]
  (fn [x]
    (map (fn [i]
           ((derivative-F (fn [xi] (f (replace-ith-vector x i xi))))
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
  (.setSensitivity tape (d+ (.sensitivity tape) sensitivity))
  (.setFanout tape (- (.fanout tape) 1))
  (when (zero? (.fanout tape))
    (let [sensitivity (.sensitivity tape)]
      (map
       (fn [factor tape] (reverse-phase! (d* sensitivity factor) tape))
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

(defn gradient-list-R [f]
  (fn [x]
    (first (second (reverse-mode
                    (fn [f x] (map f x))
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
