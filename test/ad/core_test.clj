(ns ad.core-test
  (:refer-clojure :exclude [* + - / == < > <= >= zero? number? pos? neg?])
  (:use midje.sweet)
  (:require [ad.core :refer :all]))

(facts "Basic derivatives"

  ;; f(x) = x^2
  ;; f'(x)= 2*x
  (let [f (fn [x] (* x x))
        f' (diff f)]
    (f' 0) => 0
    (f' 2) => 4
    (f' 4) => 8)

  ;; f(x) = e^x - 1.5 - atan(x)
  ;; f'(x) = e^x - 1/(1 + x^2)
  ;; NOTE multiple arguments aren't yet supported, e.g. having to nest
  ;; subtraction below
  (let [f (fn [x] (- (- (exp x) 1.5)
                     (atan x)))
        f' (diff f)]
    (f' 0) => 0.0
    (f' 1) => 2.218281828459045
    (f' 2) => 7.18905609893065))

(facts "Derivatives with a square root"
  (let [f (fn [x] (sqrt x))
        f' (diff f)]
    (f' 0) => (throws java.lang.ArithmeticException "Divide by zero")
    (f' 1) => (roughly 1/2)
    (f' 2) => 0.35355339059327373
    (f' 4) => (roughly 1/4)))

(facts "Derivatives of trigonometric functions"
  (let [f (fn [x] (sin x))
        f' (diff f)]
    (f' 0) => 1.0
    (f' (div Math/PI 2)) => (roughly 0 0.00001)
    (f' Math/PI) => -1.0)

  (let [f (fn [x] (cos x))
        f' (diff f)]
    (f' 0) => 0.0
    (f' (div Math/PI 2)) => 1.0
    (f' Math/PI) => (roughly 0.0 0.00001)
    (f' (* 3 (div Math/PI 2))) => -1.0)

  (let [f (fn [x] (tan x))
        f' (diff f)]
    (f' 0) => 1.0
    (f' Math/PI) => 1.0
    (f' (div Math/PI 2)) => 2.6670937881135714E32)

  (fact "#'log"
    (let [f (fn [x] (log x))
          f' (diff f)]
      (f 0) => Double/NEGATIVE_INFINITY
      (f' 0) => 0))

  (fact "#'exp"
    (let [f (fn [x] (exp x))
          f' (diff f)]
      (f 0.5) => 1.6487212707001282
      (f' 0.5) => 1.6487212707001282))

  (fact "#'**"
    (let [f (fn [x] (** x 2))
          f' (diff f)]
      (f 3) => 9.0
      (f' 10) => 20.0))

  (fact "#'asin"
    (let [f (fn [x] (asin x))
          f' (diff f)]
      (f 0.5) => 0.5235987755982989
      (f' 0.5) => 1.1547005383792517))

  (fact "#'acos"
    (let [f (fn [x] (acos x))
          f' (diff f)]
      (f 0.5) => 1.0471975511965979
      (f' 0.5) => 1.1547005383792517))

  (fact "hyperbolic functions"
    (let [f (fn [x] (tanh x))
          f' (diff f)]
      (f 0.5) => 0.46211715726000974
      (f' 0.5) => 0.7864477329659274)))

(facts "Computing derivatives lets us solve easily using Newton's method"
  (let [f (fn [x] (- (- (exp x) 1.5)
                     (atan x)))
        f' (diff f)
        convergent? (fn [{:keys [steps error] :or {steps 0
                                                  error 1}}]
                      (or (< error 0.00001)
                          (> steps 20)))
        next-guess (fn [{:keys [x last steps] :or {last 0
                                                  steps 0}}]
                     (let [guess (- x (div (f x) (f' x)))]
                       {:steps (inc steps)
                        :error (- last guess)
                        :last x
                        :x guess}))]
    (first (drop-while #(not (convergent? %))
                       (iterate #(next-guess %)
                                {:x -7})))
    => {:error 1.8005508195528819E-9
        :last -14.101269772739967
        :steps 7
        :x -14.101269772739967}))

(facts "gradients"
  ;; f(x,y) = x^2 + y
  ;; \del f(x,y) = [2*x, 1]
  (let [f (fn [x] (+ (* (nth x 0) (nth x 0))
                    (nth x 1)))
        f' (gradient-vector-F f)]
    (f' [2 3]) => [4 1]
    (f' [6 9]) => [12 1])

  (let [f (fn [x] (+ (* (nth x 0) (nth x 0))
                    (sin (nth x 1))))
        f' (gradient-vector-F f)]
    (f' [2 3]) => [4 -0.9899924966004454]
    (f' [6 9]) => [12 -0.9111302618846769]))

(facts "Matrix Operations"
  ;; f(x) = Ax
  ;; f'(x) = A
  (let [f (fn [x] (mmul [[1 0] [0 1]] x))
        f' (diff f)
        f'' (diff f')]
    (f [[1 2] [3 4]]) => [[1 2] [3 4]]
    (f' [[1 2] [3 4]]) => [[1 0] [0 1]]
    (f'' [[1 2] [3 4]]) => 0))

(future-facts "Satellite Trilateration"
              (let [f (fn [x t c]
                        {:pre [(= (count x) 3)
                               (= (count t) 3)
                               (number? c)]}
                        (-> (- (distance (nth x 2) (nth x 0))
                               (distance (nth x 1) (nth x 0))
                               (* c (distance (nth t 1) (nth t 2))))
                            (** 2)
                            (apply +)))
                    f' (gradient-vector-F f)]
                (f [10 15 13] [11 12 13] 1) => 0
                f' => 0))

(future-facts "Black Scholes Model"
              ;; http://www.paullegato.com/blog/black-scholes-clojure/
              ;; http://developers.opengamma.com/blog/2013/10/28/deriva-automating-algorithmic-differentiation
              (let [N (fn [x] (div 1
                                  (+ 1 (exp (- (* -0.07056 (** x 3))
                                               (* -1.5976 x))))))
                    d1 (fn [[F K T sigma]] (div (+ (div F K) (* T (div (** sigma 2) 2)))))
                    d2 (fn [[F K T sigma]] (div (- (div F K) (* T (div (** sigma 2) 2)))))
                    call (fn [[r T F K x]] (* (exp (- (* r T)))
                                             (- (* F (N (d1 [F K T x])))
                                                (* K (N (d2 [F K T x]))))))
                    put (fn [[r T F K x]] (* (exp (- (* r T)))
                                            (- (* F (N (- (d1 [F K T x]))))
                                               (* K (N (- (d2 [F K T x])))))))
                    black-model-with-sensitivities (gradient-vector-F call)]
                ;;(def black-model-with-sensitivities (∂ (bind (black-expression true) T 0.523) F K r))
                (black-model-with-sensitivities [12.3 14.3 0.03]) => nil
                (black-model-with-sensitivities [12.3 11.0 0.03]) => nil
                (black-model-with-sensitivities [12.3 11.0 0.02]) => nil))
