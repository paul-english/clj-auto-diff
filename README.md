# clj-auto-diff

This is a library for automatic differentiation in Clojure. It's based
and initiall ported from [R6RS-AD](https://github.com/qobi/R6RS-AD).

## Usage

Include it in your project using clojars,

    [clj-auto-diff "0.1.0"]

An example where being able to compute a derivative lets us solve
easily using Newton's method. Here we're solving for `x`, where `e^x -
1.5 = tan^{-1}(x)`. We do this by finding the root of `e^x - 1.5 -
tan^{-1}(x) = 0`.

```clojure
    (let [f (fn [x] (d- (d- (dexp x) 1.5)
                       (datan x)))
          f' (derivative-F f)
          convergent? (fn [{:keys [steps error] :or {steps 0
                                                    error 1}}]
                        (or (< error 0.00001)
                            (> steps 20)))
          next-guess (fn [{:keys [x last steps] :or {last 0
                                                    steps 0}}]
                       (let [guess (- x (/ (f x)
                                           (f' x)))]
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
          :x -14.101269772739967})
```

## License

Copyright © 2014 Paul English

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
